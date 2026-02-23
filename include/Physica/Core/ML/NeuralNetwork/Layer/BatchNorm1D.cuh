/*
 * Copyright 2025-2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "LayerBase.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"

namespace Physica {
    template<Scalar T>
    class BatchNorm1D;
    /**
     * Reference:
     * [1] arXiv:1502.03167; https://doi.org/10.48550/arXiv.1502.03167
     * [2] pytorch; https://pytorch.org/docs/stable/generated/torch.nn.BatchNorm1d.html
     */
    template<Scalar T>
    class device_obj<BatchNorm1D<T>> : public device_obj<LayerBase<BatchNorm1D<T>>> {
        using This = device_obj<BatchNorm1D<T>>;
        using Base = device_obj<LayerBase<BatchNorm1D<T>>>;
        using Tv = T::ValueType;
        constexpr static Tv Epsilon = std::numeric_limits<Tv>::epsilon();
        constexpr static Tv Momentum = 0.1;
    public:
        template<Scalar U>
        using MatrixND = Base::template MatrixND<U>;
    private:
        device_obj<VectorND<T>> gamma;
        device_obj<VectorND<T>> beta;
        device_obj<VectorND<Tv>> mean;
        device_obj<VectorND<Tv>> sigma;
        Tv epsilon = Epsilon;
        Tv momentum = Momentum;
    public:
        device_obj(size_t length, Tv epsilon_ = Epsilon, Tv momentum = Momentum);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) { swap(obj); return *this; }
        /* Operations */
        template<Matrix M>
        [[nodiscard]] CoDiff<device_obj<MatrixND<T>>> forward(const M& x);
        void reverse(const This& __restrict other) const noexcept;

        void step(auto& optimizer);
        void zero_grad();

        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_xavier_uniform(Tv gain);
        template<RNG R>
        void random_xavier_normal(Tv gain);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return gamma.getLength(); }
    };

    template<Scalar T>
    device_obj<BatchNorm1D<T>>::device_obj(size_t length, Tv epsilon_, Tv momentum_)
            : gamma(length, 1), beta(length, 0), mean(length, 0), sigma(length, 1), epsilon(epsilon_), momentum(momentum_) {
        assert(length > 0);
        assert(epsilon.isPositive());
        assert(momentum.isPositive() && (momentum < Tv(1)));
    }

    template<Scalar T>
    template<Matrix M>
    auto device_obj<BatchNorm1D<T>>::forward(const M& x) -> CoDiff<device_obj<MatrixND<T>>> {
        assert(getLength() == x.getRow() && "[Error]: Dimensions do not match");
        device_obj<MatrixND<Tv>> normal(x.getRow(), x.getCol());
        auto func = [normal_ = asStruct(normal),
                     x_ = asStruct(x.values()),
                     mean_ = asStruct(mean),
                     sigma_ = asStruct(sigma),
                     epsilon = epsilon,
                     momentum = momentum] __device__() mutable {
            const int i = blockIdx.x * blockDim.x + threadIdx.x;
            const auto x = x_.getDerived().row(i);
            auto& mean = mean_.getDerived()[i];
            auto& sigma = sigma_.getDerived()[i];
            auto normal = normal_.getDerived().row(i);
            mean = (Tv(1) - momentum) * mean + momentum * x.mean();
            sigma = (Tv(1) - momentum) * sigma + momentum * (x.deviation(mean) + epsilon);
            normal = (x - mean) * reciprocal(sigma);
        };
        CUDAExecutor::launch<device_obj<VectorND<T>>::MaxThreadsPerBlock>(func, gamma.makeKernelConfig());

        auto result = device_obj<MatrixND<T>>(hadamard(normal, gamma.values()) + beta.values());
        if constexpr (ReverseDiff<T>) {
            auto& result_ = co_yield std::move(result);
            if constexpr (Diffable<M>)
                x.reverse(hadamard(result_.grads(), divide(gamma.values(), sigma)));
            gamma.reverse(hadamard(result_.grads(), normal).sum_cols());
            beta.reverse(result_.grads().sum_cols());
        }
        else
            co_return std::move(result);
    }

    template<Scalar T>
    void device_obj<BatchNorm1D<T>>::reverse(const This& __restrict other) const noexcept {
        assert(this != &other && "[Error]: Self reverse is invalid");
        if constexpr (ReverseDiff<T>) {
            gamma.reverse(other.gamma.grads());
            beta.reverse(other.beta.grads());
        }
    }

    template<Scalar T>
    void device_obj<BatchNorm1D<T>>::step(auto& optimizer) {
        if constexpr (ReverseDiff<T>) {
            optimizer.step(gamma);
            optimizer.step(beta);
        }
    }

    template<Scalar T>
    void device_obj<BatchNorm1D<T>>::zero_grad() {
        if constexpr (ReverseDiff<T>) {
            gamma.zero_grad();
            beta.zero_grad();
        }
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<BatchNorm1D<T>>::random_normal() {
        beta.template random_normal<R>();
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<BatchNorm1D<T>>::random_xavier_uniform(Tv gain) {
        const auto factor = gain * sqrt(Tv(3) / Tv(getLength()));
        auto& values = beta.values();
        values.template random_uniform<R>();
        values = values * (factor * Tv(2)) - factor;
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<BatchNorm1D<T>>::random_xavier_normal(Tv gain) {
        const Tv deviation = gain * sqrt(reciprocal(Tv(getLength())));
        random_normal<R>();
        beta *= deviation;
    }

    template<Scalar T>
    void device_obj<BatchNorm1D<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is ligammaely a bug");
        gamma.swap(obj.gamma);
        beta.swap(obj.beta);
        mean.swap(obj.mean);
        sigma.swap(obj.sigma);
        epsilon.swap(obj.epsilon);
        momentum.swap(obj.momentum);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<device_obj<BatchNorm1D<T>>> {
    public:
        using ScalarType = T;
    };
}
