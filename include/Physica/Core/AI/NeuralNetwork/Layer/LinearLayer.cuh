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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType, bool WithBias>
        class Traits<device_obj<LinearLayer<ScalarType, WithBias>>> : public Traits<LinearLayer<ScalarType, WithBias>> {};
    }

    template<class ScalarType, bool WithBias>
    class device_obj<LinearLayer<ScalarType, WithBias>> : public device_obj<LayerBase<LinearLayer<ScalarType, WithBias>>> {
        static_assert(!Utils::is_device_obj<ScalarType>::value, "[Error]: Nested device_obj<> is not allowed");
        using host_obj = LinearLayer<ScalarType, WithBias>;
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
        [[nodiscard]] VectorType forward(const VectorType& x) const;
        [[nodiscard]] device_obj copy() const;

        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen);
        template<class RandomGenerator>
        void random_xavier_uniform(PlainScalar gain, RandomGenerator& gen);
        template<class RandomGenerator>
        void random_xavier_normal(PlainScalar gain, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getInputDim() const noexcept { return weights.getColumn(); }
        [[nodiscard]] __host__ __device__ size_t getOutputDim() const noexcept { return weights.getRow(); }
    };

    template<class ScalarType, bool WithBias>
    device_obj<LinearLayer<ScalarType, WithBias>>::device_obj(size_t inputDim, size_t outputDim)
            : weights(inputDim, outputDim) {
        if constexpr (WithBias)
            bias = BiasType(outputDim);
    }

    template<class ScalarType, bool WithBias>
    typename device_obj<LinearLayer<ScalarType, WithBias>>::VectorType
    device_obj<LinearLayer<ScalarType, WithBias>>::forward(const VectorType& x) const {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        if constexpr (WithBias)
            return weights * x + bias;
        else
            return weights * x;
    }

    template<class ScalarType, bool WithBias>
    device_obj<LinearLayer<ScalarType, WithBias>> device_obj<LinearLayer<ScalarType, WithBias>>::copy() const {
        This result{};
        result.weights = weights.copy();
        if constexpr (WithBias)
            result.bias = bias.copy();
        return result;
    }

    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    void device_obj<LinearLayer<ScalarType, WithBias>>::random_normal(RandomGenerator& gen) {
        weights.random_normal(gen);
        if constexpr (WithBias)
            bias.random_normal(gen);
    }
    
    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    void device_obj<LinearLayer<ScalarType, WithBias>>::random_xavier_uniform(PlainScalar gain, RandomGenerator& gen) {
        using TrivialType = typename ScalarType::TrivialType;
        const auto factor = (gain * sqrt(PlainScalar(6) / PlainScalar(getInputDim() + getOutputDim()))).getTrivial();
        std::uniform_real_distribution<TrivialType> dist(-factor, factor);
        weights.random_any(dist, gen);
        if constexpr (WithBias)
            bias.random_any(dist, gen);
    }
    
    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    void device_obj<LinearLayer<ScalarType, WithBias>>::random_xavier_normal(PlainScalar gain, RandomGenerator& gen) {
        using TrivialType = typename ScalarType::TrivialType;
        const auto deviation = (gain * sqrt(PlainScalar(2) / PlainScalar(getInputDim() + getOutputDim()))).getTrivial();
        std::normal_distribution<TrivialType> dist(0, deviation);
        weights.random_any(dist, gen);
        if constexpr (WithBias)
            bias.random_any(dist, gen);
    }
    
    template<class ScalarType, bool WithBias>
    template<class Distribution, class RandomGenerator>
    void device_obj<LinearLayer<ScalarType, WithBias>>::random_any(Distribution& dist, RandomGenerator& gen) {
        weights.random_any(dist, gen);
        if constexpr (WithBias)
            bias.random_any(dist, gen);
    }

    template<class ScalarType, bool WithBias>
    void device_obj<LinearLayer<ScalarType, WithBias>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        if constexpr (WithBias)
            bias.swap(obj.bias);
    }
}
