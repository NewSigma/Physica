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

#include "LinearCoupler.h"
#include "Physica/Core/ML/NeuralNetwork/Layer/LinearLayer.cuh"

namespace Physica {
    /**
     * Reference:
     * [1] ACM Transactions on Graphics, 38(5) 1-19 (2019); https://doi.org/10.1145/3341156
     */
    template<Scalar T>
    class device_obj<LinearCoupler<T>> {
        constexpr static int Option = MatrixOption::Row | MatrixOption::Element;
        using This = device_obj<LinearCoupler<T>>;
        using Tv = T::ValueType;
    public:
        template<Scalar U>
        using MatrixND = LinearLayer<T>::template MatrixND<U>;
    private:
        int dim;
        int numBin;
    public:
        device_obj() = default;
        device_obj(int dim_, int numBin_);
        template<Scalar U>
        device_obj(const device_obj<LinearCoupler<U>>& other);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This& obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<DNN Net>
        [[nodiscard]] CoDiff<device_obj<VectorND<T>>> forward(const Net& nn, device_obj<MatrixND<Tv>>& x) const;
        [[nodiscard]] CoDiff<device_obj<VectorND<T>>> transform(const device_obj<MatrixND<T>>& weights, device_obj<MatrixND<Tv>>& x2, int numActiveDim) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ int getDim() const noexcept { return dim; }
        [[nodiscard]] __host__ __device__ int getNumBin() const noexcept { return numBin; }
    };

    template<Scalar T>
    device_obj<LinearCoupler<T>>::device_obj(int dim_, int numBin_) : dim(dim_), numBin(numBin_) {
        assert(dim > 0);
        assert(numBin > 0);
    }

    template<Scalar T>
    template<Scalar U>
    device_obj<LinearCoupler<T>>::device_obj(const device_obj<LinearCoupler<U>>& other)
            : device_obj(other.getDim(), other.getNumBin()) {}

    template<Scalar T>
    template<DNN Net>
    auto device_obj<LinearCoupler<T>>::forward(const Net& nn, device_obj<MatrixND<Tv>>& x) const -> CoDiff<device_obj<VectorND<T>>> {
        assert(getDim() == x.getRow() && "[Error]: Dimensions do not match");
        const int dimA = getDim() / 2;
        const int dimB = dim - dimA;
        device_obj<MatrixND<Tv>> xA(getDim(), x.getCol());
        xA.zeros();
        xA.topRows(dimA) = x.topRows(dimA);
        device_obj<MatrixND<Tv>> xB = x - xA;

        auto w1 = nn.forward(xA);
        auto w2 = nn.forward(xB);
        auto lnJs = transform(w1, xB, dimB) + transform(w2, xA, dimA);
        x = xA + xB;
        if constexpr (ReverseDiff<T>) {
            auto result = co_yield lnJs.values();
            lnJs.reverse(result.grads());
        }
        else
            co_return std::move(lnJs);
    }
    /**
     * We have to make it public because of limitation of CUDA
     */
    template<Scalar T>
    auto device_obj<LinearCoupler<T>>::transform(const device_obj<MatrixND<T>>& weights, device_obj<MatrixND<Tv>>& z, int numActiveDim) const -> CoDiff<device_obj<VectorND<T>>> {
        const int numSample = z.getCol();
        auto indexes = device_obj<Array2D<size_t, Option>>(getDim(), numSample);
        auto deltas = device_obj<DenseMatrix<T, MatrixOption::Row | MatrixOption::Element>>(numSample, getDim(), 1);

        const auto factor = Tv(numBin) * Tv(1 - std::numeric_limits<Tv>::epsilon());
        auto config = device_obj<VectorND<T>>::makeKernelConfig(numSample);
        config.blocks.y = dim;
        CUDAExecutor::launch([*this,
                              weights_ = asStruct(weights),
                              z_ = asStruct(z),
                              indexes_ = asStruct(indexes),
                              deltas_ = asStruct(deltas),
                              factor] __device__() mutable {
            const int sample = blockIdx.x * blockDim.x + threadIdx.x;
            const int dim = blockIdx.y;
            const auto weights = weights_.getDerived().col(sample);
            auto z = z_.getDerived().col(sample);
            auto& indexes = indexes_.getDerived();
            auto& deltas = deltas_.getDerived();
            if (z[dim].isZero())
                return;

            const Tv tmp = z[dim] * factor;
            const int index = tmp.toMachine();
            indexes(dim, sample) = index;

            const auto grid = weights.reshape_row(getDim(), numBin);
            const auto row = grid.row(dim).values();
            const Tv lnsumexp = row.lnSumExp();
            const Tv p = softmax(row).calc(index, lnsumexp);
            deltas(sample, dim) = std::max(p, Tv(std::numeric_limits<Tv>::min()));

            Tv zi = tmp.mod() * p;
            for (int j = 0; j < index; ++j)
                zi += softmax(row).calc(j, lnsumexp);
            z[dim] = std::min(Tv(1), zi);
        }, config);

        auto expr = ln_elem(deltas);
        CoDiff<device_obj<VectorND<T>>> lnJs = expr.sum_cols() + Tv(numActiveDim) * ln(factor);
        if constexpr (ReverseDiff<T>) {
            auto tmp = co_yield lnJs.values();
            lnJs.reverse_final(tmp.grads());
            CUDAExecutor::launch([*this,
                                  weights_ = asStruct(weights),
                                  z_ = asStruct(z),
                                  indexes_ = asStruct(indexes),
                                  deltas_ = asStruct(deltas)] __device__() mutable {
                const int sample = blockIdx.x * blockDim.x + threadIdx.x;
                const int dim = blockIdx.y;
                const auto weights = weights_.getDerived().col(sample);
                const auto& indexes = indexes_.getDerived();
                const auto& deltas = deltas_.getDerived();
                auto z = z_.getDerived().col(sample);
                if (z[dim].isZero())
                    return;

                const auto grid = weights.reshape_row(getDim(), numBin);
                const int index = indexes(dim, sample);
                const Tv s = grid.row(dim).values().softmax(index);
                const Tv g = (s - square(s)) * deltas(sample, dim).grad();
                grid.row(dim).reverse(-g / Tv(numBin));
                grid(dim, index).reverse(g);
            }, config);
        }
        else
            co_return std::move(lnJs);
    }

    template<Scalar T>
    void device_obj<LinearCoupler<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(dim, obj.dim);
        std::swap(numBin, obj.numBin);
    }
}
