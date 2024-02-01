/*
 * Copyright 2023-2024 WeiBo He.
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

#include "LayerBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<class ScalarType, bool WithBias = true> class LinearLayer;

    namespace Internal {
        template<class T, bool B>
        class Traits<LinearLayer<T, B>> {
        public:
            using ScalarType = T;
            constexpr static bool WithBias = B;
            using PlainScalar = typename ScalarType::PlainScalar;
            using InputType = Vector<ScalarType>;
            using OutputType = InputType;
            constexpr static bool IsTrainMode = ScalarType::isDifferentiable;
        };
    }

    template<class ScalarType, bool WithBias>
    class LinearLayer : public LayerBase<LinearLayer<ScalarType, WithBias>> {
        using This = LinearLayer<ScalarType, WithBias>;
        using Base = LayerBase<This>;
        using typename Base::PlainScalar;
        using typename Base::InputType;
        using typename Base::OutputType;
        constexpr static int Option = MatrixOption::Row | MatrixOption::Vector;
        using MatrixType = DenseMatrix<ScalarType, Option>;
        using BiasType = typename std::conditional<WithBias, InputType, PlainStruct<void>>::type;

        MatrixType weights;
        BiasType bias;
    public:
        LinearLayer() = default;
        LinearLayer(size_t inputDim, size_t outputDim);
        template<class OtherScalar>
        LinearLayer(const LinearLayer<OtherScalar, WithBias>& layer);
        LinearLayer(const LinearLayer&) = default;
        LinearLayer(LinearLayer&&) noexcept = default;
        ~LinearLayer() = default;
        /* Operators */
        LinearLayer& operator=(LinearLayer obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] inline OutputType forward(const InputType& x) const;
        [[nodiscard]] LinearLayer copy() const;

        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen);
        template<class RandomGenerator>
        void random_xavier_uniform(PlainScalar gain, RandomGenerator& gen);
        template<class RandomGenerator>
        void random_xavier_normal(PlainScalar gain, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(LinearLayer& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getInputDim() const noexcept { return weights.getColumn(); }
        [[nodiscard]] size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const MatrixType& getWeights() const noexcept { return weights; }
        [[nodiscard]] const BiasType& getBias() const noexcept { return bias; }
    private:
        friend class device_obj<This>;
    };

    template<class ScalarType, bool WithBias>
    LinearLayer<ScalarType, WithBias>::LinearLayer(size_t inputDim, size_t outputDim)
            : weights(outputDim, inputDim) {
        if constexpr (WithBias)
            bias.resize(outputDim);
    }

    template<class ScalarType, bool WithBias>
    template<class OtherScalar>
    LinearLayer<ScalarType, WithBias>::LinearLayer(const LinearLayer<OtherScalar, WithBias>& layer)
            : weights(layer.getWeights()), bias(layer.getBias()) {}

    template<class ScalarType, bool WithBias>
    inline typename LinearLayer<ScalarType, WithBias>::OutputType LinearLayer<ScalarType, WithBias>::forward(const InputType& x) const {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        if constexpr (WithBias)
            return weights * x + bias;
        else
            return weights * x;
    }

    template<class ScalarType, bool WithBias>
    LinearLayer<ScalarType, WithBias> LinearLayer<ScalarType, WithBias>::copy() const {
        LinearLayer result{};
        result.weights = weights.copy();
        if constexpr (WithBias)
            result.bias = bias.copy();
        return result;
    }

    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    void LinearLayer<ScalarType, WithBias>::random_normal(RandomGenerator& gen) {
        weights.random_normal(gen);
        if constexpr (WithBias)
            bias.random_normal(gen);
    }

    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    void LinearLayer<ScalarType, WithBias>::random_xavier_uniform(PlainScalar gain, RandomGenerator& gen) {
        using TrivialType = typename ScalarType::TrivialType;
        const auto factor = (gain * sqrt(PlainScalar(6) / PlainScalar(getInputDim() + getOutputDim()))).getTrivial();
        std::uniform_real_distribution<TrivialType> dist(-factor, factor);
        weights.random_any(dist, gen);
        if constexpr (WithBias)
            bias.random_any(dist, gen);
    }

    template<class ScalarType, bool WithBias>
    template<class RandomGenerator>
    void LinearLayer<ScalarType, WithBias>::random_xavier_normal(PlainScalar gain, RandomGenerator& gen) {
        using TrivialType = typename ScalarType::TrivialType;
        const auto deviation = (gain * sqrt(PlainScalar(2) / PlainScalar(getInputDim() + getOutputDim()))).getTrivial();
        std::normal_distribution<TrivialType> dist(0, deviation);
        weights.random_any(dist, gen);
        if constexpr (WithBias)
            bias.random_any(dist, gen);
    }

    template<class ScalarType, bool WithBias>
    template<class Distribution, class RandomGenerator>
    void LinearLayer<ScalarType, WithBias>::random_any(Distribution& dist, RandomGenerator& gen) {
        weights.random_any(dist, gen);
        if constexpr (WithBias)
            bias.random_any(dist, gen);
    }

    template<class ScalarType, bool WithBias>
    void LinearLayer<ScalarType, WithBias>::swap(LinearLayer& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        if constexpr (WithBias)
            bias.swap(obj.bias);
    }
}
