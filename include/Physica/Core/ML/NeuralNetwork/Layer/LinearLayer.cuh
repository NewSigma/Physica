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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.cuh"

namespace Physica {
    template<Scalar T, bool WithBias>
    class device_obj<LinearLayer<T, WithBias>> : public device_obj<LayerBase<LinearLayer<T, WithBias>>> {
        static_assert(!is_device_obj<T>::value, "[Error]: Nested device_obj<> is not allowed");
        using host_obj = LinearLayer<T, WithBias>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LayerBase<host_obj>>;
        using MatrixType = host_obj::MatrixType;
        using MachineType = T::MachineType;
    public:
        using Base::IsTrain;
        using typename Base::ValueType;
        using typename Base::OutputType;
        using DiffMatrix = Diff<DenseMatrix<ValueType, host_obj::Option>, DiffMode::Reverse, T::Order>;
        using DeviceMatrix = device_obj<typename std::conditional<IsTrain, DiffMatrix, MatrixType>::type>;
        using BiasType = std::conditional<WithBias, OutputType, PlainStruct<void>>::type;
    private:
        DeviceMatrix weights;
        [[no_unique_address]] BiasType bias;
    public:
        device_obj() = default;
        template<Scalar U>
        device_obj(const device_obj<LinearLayer<U, WithBias>>& layer);
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Vector V>
        [[nodiscard]] OutputType forward(const V& x) const;

        template<RandomGenerator R>
        void random_normal();
        template<RandomGenerator R>
        void random_xavier_uniform(ValueType gain);
        template<RandomGenerator R>
        void random_xavier_normal(ValueType gain);
        template<class Distribution, RandomGenerator R>
        void random_any(Distribution& dist);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getInputDim() const noexcept { return weights.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const DeviceMatrix& getWeights() const noexcept { return weights; }
        [[nodiscard]] const BiasType& getBias() const noexcept { return bias; }
        /* Static members */
        template<RandomGenerator R>
        [[nodiscard]] static This random_normal(size_t inputDim, size_t outputDim);
        template<RandomGenerator R>
        [[nodiscard]] static This random_xavier_uniform(size_t inputDim, size_t outputDim, ValueType gain);
        template<RandomGenerator R>
        [[nodiscard]] static This random_xavier_normal(size_t inputDim, size_t outputDim, ValueType gain);
    };

    template<Scalar T, bool WithBias>
    template<Scalar U>
    device_obj<LinearLayer<T, WithBias>>::device_obj(const device_obj<LinearLayer<U, WithBias>>& layer)
            : weights(layer.getWeights())
            , bias(layer.getBias()) {}

    template<Scalar T, bool WithBias>
    template<Vector V>
    auto device_obj<LinearLayer<T, WithBias>>::forward(const V& x) const -> OutputType {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        if constexpr (WithBias)
            return weights * x + bias;
        else
            return weights * x;
    }

    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    void device_obj<LinearLayer<T, WithBias>>::random_normal() {
        weights.template random_normal<R>();
        if constexpr (WithBias)
            bias.template random_normal<R>();
    }
    
    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    void device_obj<LinearLayer<T, WithBias>>::random_xavier_uniform(ValueType gain) {
        using MachineType = T::MachineType;
        const auto factor = (gain * sqrt(ValueType(6) / ValueType(getInputDim() + getOutputDim()))).toMachine();
        std::uniform_real_distribution<MachineType> dist(-factor, factor);
        weights.template random_any<Distribution, R>(dist);
        if constexpr (WithBias)
            bias.template random_any<Distribution, R>(dist);
    }
    
    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    void device_obj<LinearLayer<T, WithBias>>::random_xavier_normal(ValueType gain) {
        const auto deviation = (gain * sqrt(ValueType(2) / ValueType(getInputDim() + getOutputDim()))).toMachine();
        std::normal_distribution<MachineType> dist(0, deviation);
        weights.template random_any<Distribution, R>(dist);
        if constexpr (WithBias)
            bias.template random_any<Distribution, R>(dist);
    }
    
    template<Scalar T, bool WithBias>
    template<class Distribution, RandomGenerator R>
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
    template<RandomGenerator R>
    device_obj<LinearLayer<T, WithBias>> device_obj<LinearLayer<T, WithBias>>::random_normal(
            size_t inputDim, size_t outputDim) {
        This result{};
        result.weights = DeviceMatrix::template random_normal<R>(outputDim, inputDim);
        if constexpr (WithBias)
            result.bias = BiasType::template random_normal<R>(outputDim);
        return result;
    }
    
    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    device_obj<LinearLayer<T, WithBias>> device_obj<LinearLayer<T, WithBias>>::random_xavier_uniform(
            size_t inputDim, size_t outputDim, ValueType gain) {
        This result{};
        const auto factor = (gain * sqrt(ValueType(6) / ValueType(inputDim + outputDim))).toMachine();
        std::uniform_real_distribution<MachineType> dist(-factor, factor);
        result.weights = DeviceMatrix::random_any(outputDim, inputDim, dist);
        if constexpr (WithBias)
            result.bias = BiasType::random_any(outputDim, dist);
        return result;
    }
    
    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    device_obj<LinearLayer<T, WithBias>> device_obj<LinearLayer<T, WithBias>>::random_xavier_normal(
            size_t inputDim, size_t outputDim, ValueType gain) {
        This result{};
        const auto deviation = (gain * sqrt(ValueType(2) / ValueType(inputDim + outputDim))).toMachine();
        std::normal_distribution<MachineType> dist(0, deviation);
        result.weights = DeviceMatrix::random_any(outputDim, inputDim, dist);
        if constexpr (WithBias)
            result.bias = BiasType::random_any(outputDim, dist);
        return result;
    }
}

namespace Physica {
    template<Scalar T, bool WithBias>
    class Traits<device_obj<LinearLayer<T, WithBias>>> : public Traits<LinearLayer<T, WithBias>> {
        using Base = Traits<LinearLayer<T, WithBias>>;
        using VectorType = VectorND<typename Base::ScalarType>;
        using ValueType = Base::ScalarType::ValueType;
        using DiffVector = Diff<VectorND<ValueType>, DiffMode::Reverse, T::Order>;
        constexpr static bool IsTrain = Base::ScalarType::isDiffable;
    public:
        using OutputType = device_obj<typename std::conditional<IsTrain, DiffVector, VectorType>::type>;
    };
}
