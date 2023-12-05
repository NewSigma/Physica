/*
 * Copyright 2023 WeiBo He.
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
    template<class ScalarType> class LinearLayer;

    namespace Internal {
        template<class T>
        class Traits<LinearLayer<T>> {
        public:
            using ScalarType = T;
        };
    }

    template<class ScalarType>
    class LinearLayer : public LayerBase<LinearLayer<ScalarType>> {
        using Base = LayerBase<LinearLayer<ScalarType>>;
        using VectorType = typename Base::VectorType;
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector>;

        MatrixType weights;
        VectorType bias;
    public:
        LinearLayer() = default;
        LinearLayer(size_t inputDim, size_t outputDim);
        template<class OtherScalar>
        LinearLayer(const LinearLayer<OtherScalar>& layer);
        LinearLayer(const LinearLayer&) = default;
        LinearLayer(LinearLayer&&) noexcept = default;
        ~LinearLayer() = default;
        /* Operators */
        LinearLayer& operator=(LinearLayer obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] inline VectorType forward(const VectorType& x) const;
        template<class Optimizer>
        inline void opt_step(const Optimizer& opt);
        [[nodiscard]] LinearLayer copy() const;

        template<class RandomGenerator>
        void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(LinearLayer& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getInputDim() const noexcept { return weights.getColumn(); }
        [[nodiscard]] size_t getOutputDim() const noexcept { return weights.getRow(); }
        [[nodiscard]] const MatrixType& getWeights() const noexcept { return weights; }
        [[nodiscard]] const VectorType& getBias() const noexcept { return bias; }
    };

    template<class ScalarType>
    LinearLayer<ScalarType>::LinearLayer(size_t inputDim, size_t outputDim) : weights(outputDim, inputDim), bias(outputDim) {}

    template<class ScalarType>
    template<class OtherScalar>
    LinearLayer<ScalarType>::LinearLayer(const LinearLayer<OtherScalar>& layer)
            : weights(layer.getWeights()), bias(layer.getBias()) {}

    template<class ScalarType>
    inline typename LinearLayer<ScalarType>::VectorType LinearLayer<ScalarType>::forward(const VectorType& x) const {
        assert(x.getLength() == getInputDim() && "[Error]: Data dim and required input dim must be equal");
        return weights * x + bias;
    }

    template<class ScalarType>
    template<class Optimizer>
    inline void LinearLayer<ScalarType>::opt_step(const Optimizer& opt) {
        opt.step(weights);
        opt.step(bias);
    }

    template<class ScalarType>
    LinearLayer<ScalarType> LinearLayer<ScalarType>::copy() const {
        LinearLayer result{};
        result.weights = weights.copy();
        result.bias = bias.copy();
        return result;
    }

    template<class ScalarType>
    template<class RandomGenerator>
    void LinearLayer<ScalarType>::random_normal(RandomGenerator& gen) {
        weights.random_normal(gen);
        bias.random_normal(gen);
    }

    template<class ScalarType>
    template<class Distribution, class RandomGenerator>
    void LinearLayer<ScalarType>::random_any(Distribution& dist, RandomGenerator& gen) {
        weights.random_any(dist, gen);
        bias.random_any(dist, gen);
    }

    template<class ScalarType>
    void LinearLayer<ScalarType>::swap(LinearLayer& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        bias.swap(obj.bias);
    }
}
