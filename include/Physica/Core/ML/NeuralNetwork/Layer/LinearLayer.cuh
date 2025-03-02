/*
 * Copyright 2024-2025 Weibo He.
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

#include "LinearLayer.h"
#include "LayerBase.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.cuh"

namespace Physica {
    template<Scalar T, bool WithBias>
    class device_obj<LinearLayer<T, WithBias>> : public device_obj<LayerBase<LinearLayer<T, WithBias>>> {
        using host_obj = LinearLayer<T, WithBias>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerBase<host_obj>>;

        using HostMatrix = host_obj::MatrixType;
        using MachineType = T::MachineType;
        using MatrixType = device_obj<DenseMatrix<T, host_obj::Option>>;
        using BiasType = std::conditional<WithBias, device_obj<VectorND<T>>, PlainStruct<void>>::type;
    private:
        using Tv = T::ValueType;

        MatrixType weights;
        [[no_unique_address]] BiasType bias;
    public:
        device_obj() = default;
        device_obj(size_t inputDim, size_t outputDim);
        template<Scalar U>
        device_obj(const device_obj<LinearLayer<U, WithBias>>& layer);
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Vector V>
        [[nodiscard]] inline CoDiff<device_obj<VectorND<T>>> forward(const V& x) const;
        void reverse(const This& __restrict other) const noexcept;

        template<class Optimizer>
        void step(Optimizer& opt);
        void zero_grad();

        void resize(size_t inputDim, size_t outputDim);
        [[nodiscard]] inline host_obj toHost() const;
        [[nodiscard]] inline host_obj toHostAsync() const;
        inline void toHost(host_obj& obj) const;
        inline void toHostAsync(host_obj& obj) const;

        void fill(Tv x);
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_xavier_uniform(Tv gain);
        template<RNG R>
        void random_xavier_normal(Tv gain);
        template<RNG R, class Distribution>
        void random_any(Distribution& dist);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getInputDim() const noexcept { return weights.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const auto& getWeights() const noexcept { return weights; }
        [[nodiscard]] const auto& getBias() const noexcept requires(WithBias) { return bias; }
    private:
        template<Vector V>
        [[nodiscard]] inline CoDiff<device_obj<VectorND<T>>> forward_impl(const V& x) const requires(CUDA<V>);
    };

    template<Scalar T, bool WithBias>
    device_obj<LinearLayer<T, WithBias>>::device_obj(size_t inputDim, size_t outputDim) {
        resize(inputDim, outputDim);
    }

    template<Scalar T, bool WithBias>
    template<Scalar U>
    device_obj<LinearLayer<T, WithBias>>::device_obj(const device_obj<LinearLayer<U, WithBias>>& layer) {
        if constexpr (Diffable<T>) {
            weights = layer.getWeights();
            if constexpr (WithBias)
                bias = layer.getBias();
        }
        else {
            weights = layer.getWeights().values();
            if constexpr (WithBias)
                bias = layer.getBias().values();
        }
    }

    template<Scalar T, bool WithBias>
    template<Vector V>
    auto device_obj<LinearLayer<T, WithBias>>::forward(const V& x) const -> CoDiff<device_obj<VectorND<T>>> {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        if constexpr (CUDA<V>)
            return forward_impl(x);
        else
            return forward_impl(x.toDevice());
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::reverse(const This& __restrict other) const noexcept {
        assert(this != &other && "[Error]: Self reverse is invalid");
        assert(getInputDim() == other.getInputDim());
        assert(getOutputDim() == other.getOutputDim());
        if constexpr (ReverseDiff<T>) {
            weights.reverse(other.weights.grads());
            if constexpr (WithBias)
                bias.reverse(other.bias.grads());
        }
    }

    template<Scalar T, bool WithBias>
    template<class Optimizer>
    void device_obj<LinearLayer<T, WithBias>>::step(Optimizer& opt) {
        if constexpr (ReverseDiff<T>) {
            opt.step(weights);
            if constexpr (WithBias)
                opt.step(bias);
        }
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::zero_grad() {
        if constexpr (ReverseDiff<T>) {
            weights.grads().zeros();
            if constexpr (WithBias)
                bias.grads().zeros();
        }
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::resize(size_t inputDim, size_t outputDim) {
        weights.resize(outputDim, inputDim);
        if constexpr (WithBias)
            bias.resize(outputDim);
    }

    template<Scalar T, bool WithBias>
    inline auto device_obj<LinearLayer<T, WithBias>>::toHost() const -> host_obj {
        host_obj result = toHostAsync();
        CUDAExecutor::wait();
        return result;
    }

    template<Scalar T, bool WithBias>
    inline auto device_obj<LinearLayer<T, WithBias>>::toHostAsync() const -> host_obj {
        host_obj result(getInputDim(), getOutputDim());
        toHostAsync(result);
        return result;
    }

    template<Scalar T, bool WithBias>
    inline void device_obj<LinearLayer<T, WithBias>>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAExecutor::wait();
    }

    template<Scalar T, bool WithBias>
    inline void device_obj<LinearLayer<T, WithBias>>::toHostAsync(host_obj& obj) const {
        weights.toHostAsync(obj.getWeights());
        if constexpr (WithBias)
            bias.toHostAsync(obj.getBias());
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::fill(Tv x) {
        weights = x;
        if constexpr (WithBias)
            bias = x;
    }

    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_normal() {
        weights.template random_normal<R>();
        if constexpr (WithBias)
            bias.template random_normal<R>();
    }
    
    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_xavier_uniform(Tv gain) {
        const auto factor = (gain * sqrt(Tv(6) / Tv(getInputDim() + getOutputDim()))).toMachine();
        std::uniform_real_distribution<MachineType> dist(-factor, factor);
        weights.template random_any<R, decltype(dist)>(dist);
        if constexpr (WithBias)
            bias.template random_any<R, decltype(dist)>(dist);
    }
    
    template<Scalar T, bool WithBias>
    template<RNG R>
    void device_obj<LinearLayer<T, WithBias>>::random_xavier_normal(Tv gain) {
        const auto deviation = (gain * sqrt(Tv(2) / Tv(getInputDim() + getOutputDim()))).toMachine();
        std::normal_distribution<MachineType> dist(0, deviation);
        weights.template random_any<R, decltype(dist)>(dist);
        if constexpr (WithBias)
            bias.template random_any<R, decltype(dist)>(dist);
    }
    
    template<Scalar T, bool WithBias>
    template<RNG R, class Distribution>
    void device_obj<LinearLayer<T, WithBias>>::random_any(Distribution& dist) {
        weights.random_any(dist);
        if constexpr (WithBias)
            bias.random_any(dist);
    }

    template<Scalar T, bool WithBias>
    void device_obj<LinearLayer<T, WithBias>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        if constexpr (WithBias)
            bias.swap(obj.bias);
    }

    template<Scalar T, bool WithBias>
    template<Vector V>
    inline CoDiff<device_obj<VectorND<T>>> device_obj<LinearLayer<T, WithBias>>::forward_impl(const V& x) const requires(CUDA<V>) {
        if constexpr (ReverseDiff<T>) {
            auto expr1 = weights * x;
            if constexpr (WithBias) {
                auto expr2 = expr1 + bias;
                const auto result = co_yield expr2.values();
                expr2.reverse(result.grads());
            }
            else {
                const auto result = co_yield expr1.values();
                expr1.reverse(result.grads());
            }
        }
        else {
            if constexpr (WithBias)
                co_return weights * x + bias;
            else
                co_return weights * x;
        }
    }
}

namespace Physica {
    template<Scalar T, bool WithBias>
    class Traits<device_obj<LinearLayer<T, WithBias>>> : public Traits<LinearLayer<T, WithBias>> {};
}
