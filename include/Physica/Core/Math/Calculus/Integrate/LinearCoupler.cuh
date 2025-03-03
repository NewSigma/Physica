/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Utils/Container/Array.cuh"
#include "LinearCoupler.h"

namespace Physica {
    /**
     * Reference:
     * [1] ACM Transactions on Graphics, 38(5) 1-19 (2019); https://doi.org/10.1145/3341156
     */
    template<Scalar T>
    class device_obj<LinearCoupler<T>> {
        using This = device_obj<LinearCoupler<T>>;
        using Tv = T::ValueType;
    private:
        device_obj<VectorND<Tv>> mask;
        int numBin;
    public:
        device_obj() = default;
        device_obj(int dim, int numBin_);
        template<Scalar U>
        device_obj(const device_obj<LinearCoupler<U>>& other);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This& obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<DNN Net>
        [[nodiscard]] CoDiff<T> forward(const Net& nn, device_obj<VectorND<Tv>>& x) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getDim() const noexcept { return mask.getLength(); }
        [[nodiscard]] int getNumBin() const noexcept { return numBin; }
    private:
        CoDiff<T> transform(const device_obj<VectorND<T>>& weights, device_obj<VectorND<Tv>>& x2, int numActiveDim) const;
    };

    template<Scalar T>
    device_obj<LinearCoupler<T>>::device_obj(int dim, int numBin_) : mask(dim), numBin(numBin_) {
        assert(dim > 0);
        assert(numBin > 0);
        mask.zeros();
        mask.head(dim / 2) = Tv(1);
    }

    template<Scalar T>
    template<Scalar U>
    device_obj<LinearCoupler<T>>::device_obj(const device_obj<LinearCoupler<U>>& other)
            : device_obj(other.getDim(), other.getNumBin()) {}

    template<Scalar T>
    template<DNN Net>
    auto device_obj<LinearCoupler<T>>::forward(const Net& nn, device_obj<VectorND<Tv>>& x) const -> CoDiff<T> {
        assert(getDim() == x.getLength() && "[Error]: Dimensions do not match");
        device_obj<VectorND<Tv>> xA = hadamard(x, mask);
        device_obj<VectorND<Tv>> xB = x - xA;
        auto w1 = nn.forward(xA);
        auto w2 = nn.forward(xB);
        auto lnJ = transform(w1, xB) + transform(w2, xA);
        x = xA + xB;
        if constexpr (ReverseDiff<T>) {
            auto result = co_yield lnJ.value();
            lnJ.reverse(result.grad());
        }
        else
            co_return std::move(lnJ);
    }

    template<Scalar T>
    void device_obj<LinearCoupler<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mask.swap(obj.mask);
        std::swap(numBin, obj.numBin);
    }

    template<Scalar T>
    auto device_obj<LinearCoupler<T>>::transform(const device_obj<VectorND<T>>& weights, device_obj<VectorND<Tv>>& z, int numActiveDim) const -> CoDiff<T> {
        const size_t length = z.getLength();
        device_obj<Array<size_t>> indexes(length);
        device_obj<VectorND<T>> deltas(length, 1);

        const auto factor = Tv(numBin) * Tv(1 - std::numeric_limits<Tv>::epsilon());
        const auto config = z.makeKernelConfig();
        CUDAExecutor::launch([*this, weights_ = asStruct(weights), z_ = asStruct(z), factor, numActiveDim = Tv(numActiveDim)] __device__() mutable {
            const int i = blockIdx.x * blockDim.x + threadIdx.x;
            const auto& weights = weights_.getDerived();
            auto& z = z_.getDerived();
            
            VectorND<Tv> prob(numBin);

            if (z[i].isZero())
                return;
            const auto grid = weights.reshape_col(numBin, getDim());
            prob = softmax(grid.col(i).values());

            const Tv tmp = z[i] * factor;
            const size_t index = tmp.toMachine();
            indexes[i] = index;
            deltas[i] = std::max(prob[index], Tv(std::numeric_limits<Tv>::min()));
            Tv zi = tmp.mod() * prob[index];
            if (index > 0)
                zi = std::min(Tv(1), zi + prob.head(index).sum());
            z[i] = zi;
        }, config.first, config.second);

        auto expr = ln(deltas);
        auto lnJ = expr.sum() + numActiveDim * ln(factor);
        if constexpr (ReverseDiff<T>) {
            auto tmp = co_yield lnJ.value();
            lnJ.reverse_final(tmp.grad());
            const auto config = z.makeKernelConfig();
            CUDAExecutor::launch([*this, weights_ = asStruct(weights), z_ = asStruct(z)] __device__() mutable {
                const int i = blockIdx.x * blockDim.x + threadIdx.x;
                const auto& weights = weights_.getDerived();
                auto& z = z_.getDerived();
                if (z[i].isZero())
                    return;

                const auto grid = weights.reshape_col(numBin, getDim());
                Tv s = grid.col(i).values().softmax(indexes[i]);
                Tv g = (s - square(s)) * deltas[i].grad();
                grid.col(i).reverse(-g / Tv(numBin));
                grid(indexes[i], i).reverse(g);
            }, config.first, config.second);
        }
        else
            co_return std::move(lnJ);
    }
}
