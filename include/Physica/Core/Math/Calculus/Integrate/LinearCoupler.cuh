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
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<DNN Net>
        [[nodiscard]] CoDiff<device_obj<VectorND<T>>> forward(Net& nn, device_obj<MatrixND<Tv>>& x) const;
        [[nodiscard]] CoDiff<device_obj<VectorND<T>>> transform(const device_obj<MatrixND<T>>& weights, device_obj<MatrixND<Tv>>& x2) const;

        template<RNG R>
        void random_shuffle();

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ int getDim() const noexcept { return mask.getLength(); }
        [[nodiscard]] __host__ __device__ int getNumBin() const noexcept { return numBin; }
    private:
        Tv calcFactor() const noexcept { return Tv(numBin) * Tv(1 - std::numeric_limits<Tv>::epsilon()); }
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
    auto device_obj<LinearCoupler<T>>::forward(Net& nn, device_obj<MatrixND<Tv>>& x) const -> CoDiff<device_obj<VectorND<T>>> {
        assert(getDim() == x.getRow() && "[Error]: Dimensions do not match");
        device_obj<MatrixND<Tv>> xA = hadamard(x, mask);
        device_obj<MatrixND<Tv>> xB = x - xA;

        auto w1 = nn.forward(xA);
        auto w2 = nn.forward(xB);
        auto lnJs1 = transform(w1, xB);
        auto lnJs2 = transform(w2, xA);
        x = xA + xB;
        if constexpr (ReverseDiff<T>) {
            auto result = co_yield lnJs1.values() + lnJs2.values() + Tv(getDim()) * ln(calcFactor());
            lnJs1.reverse(result.grads());
            lnJs2.reverse(result.grads());
        }
        else
            co_return lnJs1 + lnJs2;
    }
    /**
     * We have to make it public because of limitation of CUDA
     */
    template<Scalar T>
    auto device_obj<LinearCoupler<T>>::transform(const device_obj<MatrixND<T>>& weights, device_obj<MatrixND<Tv>>& z) const -> CoDiff<device_obj<VectorND<T>>> {
        const int numSample = z.getCol();
        auto indices = device_obj<Array2D<size_t, Option>>(getDim(), numSample);
        auto deltas = device_obj<DenseMatrix<T, MatrixOption::Row | MatrixOption::Element>>(numSample, getDim(), 1);

        auto config = device_obj<VectorND<T>>::makeKernelConfig(numSample);
        config.blocks.y = getDim();
        auto fwd = [dim = getDim(),
                    numBin = numBin,
                    weights_ = asStruct(weights),
                    z_ = asStruct(z),
                    indices_ = asStruct(indices),
                    deltas_ = asStruct(deltas),
                    factor = calcFactor()] __device__() mutable {
            const int sample = blockIdx.x * blockDim.x + threadIdx.x;
            const int i = blockIdx.y;
            const auto weights = weights_.getDerived().col(sample);
            auto z = z_.getDerived().col(sample);
            auto& indices = indices_.getDerived();
            auto& deltas = deltas_.getDerived();
            if (z[i].isZero())
                return;

            const Tv tmp = z[i] * factor;
            const int index = tmp.toMachine();
            indices(i, sample) = index;

            const auto grid = weights.reshape_row(dim, numBin);
            const auto row = grid.row(i).values();
            const Tv lnsumexp = row.lnSumExp();
            const Tv p = softmax(row).calc(index, lnsumexp);
            deltas(sample, i) = std::max(p, Tv(std::numeric_limits<Tv>::min()));

            Tv zi = tmp.mod() * p;
            for (int j = 0; j < index; ++j)
                zi += softmax(row).calc(j, lnsumexp);
            z[i] = std::min(Tv(1), zi);
        };
        CUDAExecutor::launch<decltype(fwd), device_obj<VectorND<T>>::MaxThreadPerBlock>(fwd, config);

        auto expr = ln_elem(deltas);
        CoDiff<device_obj<VectorND<T>>> lnJs = expr.sum_cols();
        if constexpr (ReverseDiff<T>) {
            auto tmp = co_yield lnJs.values();
            lnJs.reverse_final(tmp.grads());
            auto func = [dim = getDim(),
                         numBin = numBin,
                         weights_ = asStruct(weights),
                         z_ = asStruct(z),
                         indices_ = asStruct(indices),
                         deltas_ = asStruct(deltas)] __device__() mutable {
                const int sample = blockIdx.x * blockDim.x + threadIdx.x;
                const int i = blockIdx.y;
                const auto weights = weights_.getDerived().col(sample);
                const auto& indices = indices_.getDerived();
                const auto& deltas = deltas_.getDerived();
                auto z = z_.getDerived().col(sample);
                if (z[i].isZero())
                    return;

                const auto grid = weights.reshape_row(dim, numBin);
                const int index = indices(i, sample);
                const auto row = grid.row(i);
                const Tv s = row.values().softmax(index);
                const Tv g = (s - square(s)) * deltas(sample, i).grad();
                row.reverse(-g / Tv(numBin));
                row[index].reverse(g);
            };
            CUDAExecutor::launch<decltype(func), device_obj<VectorND<T>>::MaxThreadPerBlock>(func, config);
        }
        else
            co_return std::move(lnJs);
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<LinearCoupler<T>>::random_shuffle() {
        VectorND<Tv> mask_shuffle(getDim());
        mask_shuffle.zeros();
        mask_shuffle.head(getDim() / 2) = Tv(1);
        std::shuffle(mask_shuffle.begin(), mask_shuffle.end(), R::getInstance());
        mask_shuffle.toDeviceAsync(mask);
    }

    template<Scalar T>
    void device_obj<LinearCoupler<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mask.swap(obj.mask);
        std::swap(numBin, obj.numBin);
    }
}
