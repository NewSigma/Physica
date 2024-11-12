/*
 * Copyright 2024 Weibo He.
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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.cuh>

namespace Physica::Core {
    template<class ScalarType, bool WithBias>
    class device_obj<LinearLayer<ScalarType, WithBias>> : public device_obj<LayerBase<LinearLayer<ScalarType, WithBias>>> {
        static_assert(!is_device_obj<ScalarType>::value, "[Error]: Nested device_obj<> is not allowed");
        using host_obj = LinearLayer<ScalarType, WithBias>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerBase<host_obj>>;
        using MatrixType = typename host_obj::MatrixType;
        using MachineType = typename ScalarType::MachineType;
    public:
        using Base::IsTrainMode;
        using typename Base::ValueType;
        using typename Base::InputType;
        using typename Base::OutputType;
        using DiffMatrix = Diff<DenseMatrix<ValueType, host_obj::Option>, DiffMode::Reverse, ScalarType::Order>;
        using DeviceMatrix = device_obj<typename std::conditional<IsTrainMode, DiffMatrix, MatrixType>::type>;
        using BiasType = typename std::conditional<WithBias, InputType, PlainStruct<void>>::type;
    private:
        DeviceMatrix weights;
        BiasType bias;
    public:
        device_obj() = default;
        template<class OtherScalar>
        device_obj(const device_obj<LinearLayer<OtherScalar, WithBias>>& layer);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] OutputType forward(const InputType& x) const;
        [[nodiscard]] device_obj copy() const;

        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen);
        template<class RandomGenerator>
        void random_xavier_uniform(ValueType gain, RandomGenerator& gen);
        template<class RandomGenerator>
        void random_xavier_normal(ValueType gain, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getInputDim() const noexcept { return weights.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const DeviceMatrix& getWeights() const noexcept { return weights; }
        [[nodiscard]] const BiasType& getBias() const noexcept { return bias; }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] static This random_normal(size_t inputDim, size_t outputDim, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static This random_xavier_uniform(size_t inputDim, size_t outputDim, ValueType gain, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static This random_xavier_normal(size_t inputDim, size_t outputDim, ValueType gain, RandomGenerator& gen);
    };

    template<class ScalarType, bool WithBias>
    template<class OtherScalar>
    device_obj<LinearLayer<ScalarType, WithBias>>::device_obj(const device_obj<LinearLayer<OtherScalar, WithBias>>& layer)
            : weights(layer.getWeights())
            , bias(layer.getBias()) {}

    template<class ScalarType, bool WithBias>
    typename device_obj<LinearLayer<ScalarType, WithBias>>::OutputType
    device_obj<LinearLayer<ScalarType, WithBias>>::forward(const InputType& x) const {
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
    void device_obj<LinearLayer<ScalarType, WithBias>>::random_xavier_uniform(ValueType gain, RandomGenerator& gen) {
        using MachineType = typename ScalarType::MachineType;
        const auto factor = (gain * sqrt(ValueType(6) / ValueType(getInputDim() + getOutputDim()))).toMachine();
        std::uniform_real_distribution<MachineType> dist(-factor, factor);
        weights.random_any(dist, gen);
        if constexpr (WithBias)
            bias.random_any(dist, gen);
    }
    
    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    void device_obj<LinearLayer<ScalarType, WithBias>>::random_xavier_normal(ValueType gain, RandomGenerator& gen) {
        const auto deviation = (gain * sqrt(ValueType(2) / ValueType(getInputDim() + getOutputDim()))).toMachine();
        std::normal_distribution<MachineType> dist(0, deviation);
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

    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    device_obj<LinearLayer<ScalarType, WithBias>> device_obj<LinearLayer<ScalarType, WithBias>>::random_normal(
            size_t inputDim, size_t outputDim, RandomGenerator& gen) {
        This result{};
        result.weights = DeviceMatrix::random_normal(outputDim, inputDim, gen);
        if constexpr (WithBias)
            result.bias = BiasType::random_normal(outputDim, gen);
        return result;
    }
    
    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    device_obj<LinearLayer<ScalarType, WithBias>> device_obj<LinearLayer<ScalarType, WithBias>>::random_xavier_uniform(
            size_t inputDim, size_t outputDim, ValueType gain, RandomGenerator& gen) {
        This result{};
        const auto factor = (gain * sqrt(ValueType(6) / ValueType(inputDim + outputDim))).toMachine();
        std::uniform_real_distribution<MachineType> dist(-factor, factor);
        result.weights = DeviceMatrix::random_any(outputDim, inputDim, dist, gen);
        if constexpr (WithBias)
            result.bias = BiasType::random_any(outputDim, dist, gen);
        return result;
    }
    
    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    device_obj<LinearLayer<ScalarType, WithBias>> device_obj<LinearLayer<ScalarType, WithBias>>::random_xavier_normal(
            size_t inputDim, size_t outputDim, ValueType gain, RandomGenerator& gen) {
        This result{};
        const auto deviation = (gain * sqrt(ValueType(2) / ValueType(inputDim + outputDim))).toMachine();
        std::normal_distribution<MachineType> dist(0, deviation);
        result.weights = DeviceMatrix::random_any(outputDim, inputDim, dist, gen);
        if constexpr (WithBias)
            result.bias = BiasType::random_any(outputDim, dist, gen);
        return result;
    }
}

namespace Physica {
    using namespace Core;

    template<class T, bool WithBias>
    class Traits<Core::device_obj<LinearLayer<T, WithBias>>> : public Traits<LinearLayer<T, WithBias>> {
        using Base = Traits<LinearLayer<T, WithBias>>;
        using VectorType = Vector<typename Base::ScalarType>;
        using ValueType = typename Base::ScalarType::ValueType;
        using DiffVector = Diff<Vector<ValueType>, DiffMode::Reverse, T::Order>;
        constexpr static bool IsTrainMode = Base::ScalarType::isDifferentiable;
    public:
        using InputType = Core::device_obj<typename std::conditional<IsTrainMode, DiffVector, VectorType>::type>;
        using OutputType = InputType;
    };
}
