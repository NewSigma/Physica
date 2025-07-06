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

#include "Physica/Core/ML/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.cuh"

namespace Physica {
    template<Scalar T>
    class VegasLayer;

    template<Scalar T>
    class device_obj<VegasLayer<T>> {
        using host_obj = VegasLayer<T>;
        using This = device_obj<host_obj>;
        using Tv = T::ValueType;
    public:
        template<Scalar U>
        using MatrixND = LinearLayer<T>::template MatrixND<U>;
    private:
        device_obj<VectorND<T>> weights;
        int dim;
    public:
        device_obj(int dim_, int numBin);
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Matrix M>
        [[nodiscard]] CoDiff<device_obj<MatrixND<T>>> forward(const M& x) const requires(CUDA<M>);

        void step(auto& optimizer);
        void zero_grad();

        const auto read(const H5Loc& loc, const char* name);
        auto write(H5Loc& loc, const char* name) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getWeights() const noexcept { return weights; }
        [[nodiscard]] int getDim() const noexcept { return dim; }
        [[nodiscard]] int getNumBin() const noexcept { return weights.getLength() / getDim(); }
    };

    template<Scalar T>
    device_obj<VegasLayer<T>>::device_obj(int dim_, int numBin) : weights(dim_ * numBin), dim(dim_) {
        weights.zeros();
    }

    template<Scalar T>
    template<Matrix M>
    auto device_obj<VegasLayer<T>>::forward(const M& x) const -> CoDiff<device_obj<MatrixND<T>>> requires(CUDA<M>) {
        assert(x.getRow() == weights.getLength());
        CoDiff<device_obj<MatrixND<T>>> y = x + weights;
        if constexpr (ReverseDiff<T>) {
            auto co_y = co_yield y.values();
            y.reverse(co_y.grads());
        }
        else
            co_return std::move(y);
    }

    template<Scalar T>
    void device_obj<VegasLayer<T>>::step(auto& optimizer) {
        optimizer.step(weights);

        auto func = [v_ = asStruct(weights.values()), numDim = dim] __device__() mutable {
            extern __shared__ Tv buffer[];
            const int dim = blockIdx.x;
            auto& v = v_.getDerived();
            auto mat = v.reshape_row(numDim, v.getLength() / numDim);
            auto row = mat.row(dim);
            const Tv mean = row.mean(threadIdx.x, blockDim.x, buffer);
            const int length = row.getLength();
            for (int i = threadIdx.x; i < length; i += blockDim.x)
                row[i] -= mean;
        };
        const int numThreads = std::min(getNumBin(), weights.MaxThreadPerBlock);
        CUDAExecutor::launch(func, KernelConfig(getDim(), numThreads), numThreads * sizeof(Tv));
    }

    template<Scalar T>
    void device_obj<VegasLayer<T>>::zero_grad() {
        weights.zero_grad();
    }

#ifdef PHYSICA_HDF5
template<Scalar T>
    const auto device_obj<VegasLayer<T>>::read(const H5Loc& loc, const char* name) {
        return weights.values().read(loc, name);
    }

    template<Scalar T>
    auto device_obj<VegasLayer<T>>::write(H5Loc& loc, const char* name) const {
        return weights.values().write(loc, name);
    }
#endif

    template<Scalar T>
    void device_obj<VegasLayer<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        std::swap(dim, obj.dim);
    }
}
