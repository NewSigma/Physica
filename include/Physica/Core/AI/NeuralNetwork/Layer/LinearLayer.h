/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/PlainStruct.h"
#include "LayerBase.h"

namespace Physica::Core {
    template<Scalar T, bool WithBias = true>
    class LinearLayer : public LayerBase<LinearLayer<T, WithBias>> {
        using This = LinearLayer<T, WithBias>;
        using Base = LayerBase<This>;
        using typename Base::InputType;
        using typename Base::OutputType;
        using typename Base::ValueType;
        constexpr static int Option = MatrixOption::Row | MatrixOption::Vector;
        using MatrixType = DenseMatrix<T, Option>;
        using BiasType = std::conditional<WithBias, InputType, PlainStruct<void>>::type;
    public:
        using device_obj_type = device_obj<This>;
    private:
        MatrixType weights;
        [[no_unique_address]] BiasType bias;
    public:
        LinearLayer() = default;
        LinearLayer(size_t inputDim, size_t outputDim);
        template<Scalar U>
        LinearLayer(const LinearLayer<U, WithBias>& layer);
        LinearLayer(const This&) = default;
        LinearLayer(This&&) noexcept = default;
        ~LinearLayer() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] inline OutputType forward(const InputType& x) const;
        [[nodiscard]] LinearLayer copy() const;

        template<RandomGenerator R>
        void random_normal();
        template<RandomGenerator R>
        void random_xavier_uniform(ValueType gain);
        template<RandomGenerator R>
        void random_xavier_normal(ValueType gain);
        template<class Distribution, RandomGenerator R>
        void random_any(Distribution& dist);
        void swap(LinearLayer& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getInputDim() const noexcept { return weights.getCol(); }
        [[nodiscard]] size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const MatrixType& getWeights() const noexcept { return weights; }
        [[nodiscard]] const BiasType& getBias() const noexcept { return bias; }
    private:
        friend class device_obj<This>;
    };

    template<Scalar T, bool WithBias>
    LinearLayer<T, WithBias>::LinearLayer(size_t inputDim, size_t outputDim)
            : weights(outputDim, inputDim) {
        if constexpr (WithBias)
            bias.resize(outputDim);
    }

    template<Scalar T, bool WithBias>
    template<Scalar U>
    LinearLayer<T, WithBias>::LinearLayer(const LinearLayer<U, WithBias>& layer)
            : weights(layer.getWeights()), bias(layer.getBias()) {}

    template<Scalar T, bool WithBias>
    inline LinearLayer<T, WithBias>::OutputType LinearLayer<T, WithBias>::forward(const InputType& x) const {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        if constexpr (WithBias)
            return weights * x + bias;
        else
            return weights * x;
    }

    template<Scalar T, bool WithBias>
    LinearLayer<T, WithBias> LinearLayer<T, WithBias>::copy() const {
        LinearLayer result{};
        result.weights = weights.copy();
        if constexpr (WithBias)
            result.bias = bias.copy();
        return result;
    }

    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    void LinearLayer<T, WithBias>::random_normal() {
        weights.template random_normal<R>();
        if constexpr (WithBias)
            bias.template random_normal<R>();
    }

    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    void LinearLayer<T, WithBias>::random_xavier_uniform(ValueType gain) {
        using MachineType = T::MachineType;
        const auto factor = (gain * sqrt(ValueType(6) / ValueType(getInputDim() + getOutputDim()))).toMachine();
        std::uniform_real_distribution<MachineType> dist(-factor, factor);
        weights.template random_any<R>(dist);
        if constexpr (WithBias)
            bias.template random_any<R>(dist);
    }

    template<Scalar T, bool WithBias>
    template<RandomGenerator R>
    void LinearLayer<T, WithBias>::random_xavier_normal(ValueType gain) {
        using MachineType = T::MachineType;
        const auto deviation = (gain * sqrt(ValueType(2) / ValueType(getInputDim() + getOutputDim()))).toMachine();
        std::normal_distribution<MachineType> dist(0, deviation);
        weights.template random_any<decltype(dist), R>(dist);
        if constexpr (WithBias)
            bias.template random_any<decltype(dist), R>(dist);
    }

    template<Scalar T, bool WithBias>
    template<class Distribution, RandomGenerator R>
    void LinearLayer<T, WithBias>::random_any(Distribution& dist) {
        weights.template random_any<Distribution, R>(dist);
        if constexpr (WithBias)
            bias.template random_any<Distribution, R>(dist);
    }

    template<Scalar T, bool WithBias>
    void LinearLayer<T, WithBias>::swap(LinearLayer& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        if constexpr (WithBias)
            bias.swap(obj.bias);
    }
}

namespace Physica {
    template<Scalar T, bool B>
    class Traits<Core::LinearLayer<T, B>> {
    public:
        using ScalarType = T;
        constexpr static bool WithBias = B;

        using InputType = Core::VectorND<T>;
        using OutputType = InputType;
    };
}
