/*
 * Copyright 2024 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.cuh"

namespace Physica::Core {
    template<class ScalarType, bool WithBias>
    class device_obj<Linear<ScalarType, WithBias>> : device_obj<LayerBase<Linear<ScalarType, WithBias>>> {
        using host_obj = Linear<ScalarType, WithBias>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerBase<host_obj>>;
        using MatrixType = typename host_obj::MatrixType;
    public:
        using Base::IsTrainMode;
        using typename Base::VectorType;
        using typename Base::PlainScalar;
        using DiffMatrix = Differentiable<DenseMatrix<PlainScalar, host_obj::Option>, DiffMode::Reverse>;
        using DeviceMatrix = device_obj<typename std::conditional<IsTrainMode, DiffMatrix, MatrixType>::type>;
        using BiasType = typename std::conditional<WithBias, VectorType, PlainStruct<void>>::type;
    private:
        DeviceMatrix weights;
        BiasType bias;
    public:
        device_obj() = default;
        device_obj(size_t inputDim, size_t outputDim);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] VectorType forward(const VectorType& input) const;
        [[nodiscard]] device_obj copy() const;
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getInputDim() const noexcept { return weights.getColumn(); }
        [[nodiscard]] __host__ __device__ size_t getOutputDim() const noexcept { return weights.getRow(); }
    }

    template<class ScalarType, bool WithBias>
    device_obj<Linear<ScalarType, WithBias>>::device_obj(size_t inputDim, size_t outputDim)
            : weights(inputDim, outputDim) {
        if constexpr (WithBias)
            bias = BiasType(outputDim);
    }

    template<class ScalarType, bool WithBias>
    typename device_obj<Linear<ScalarType, WithBias>>::VectorType
    device_obj<Linear<ScalarType, WithBias>>::forward(const VectorType& input) const {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        if constexpr (WithBias)
            return weights * x + bias;
        else
            return weights * x;
    }

    template<class ScalarType, bool WithBias>
    device_obj<Linear<ScalarType, WithBias>> device_obj<Linear<ScalarType, WithBias>>::copy() const {
        This result{};
        result.weights = weights.copy();
        if constexpr (WithBias)
            result.bias = bias.copy();
        return result;
    }

    template<class ScalarType, bool WithBias>
    void device_obj<Linear<ScalarType, WithBias>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        if constexpr (WithBias)
            bias.swap(obj.bias);
    }
}
