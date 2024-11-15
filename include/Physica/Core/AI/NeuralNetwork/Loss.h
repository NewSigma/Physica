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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h>

namespace Physica::Core {
    template<class ScalarType>
    class Loss {
        static_assert(!is_device_obj<ScalarType>::value, "[Error]: Include corresponding *.cuh file to enable CUDA support");
    public:
        constexpr static bool IsTrainMode = ScalarType::isDifferentiable;
        using ValueType = typename ScalarType::ValueType;
        using LossType = ScalarType;
    private:
        using VectorType = VectorND<ScalarType>;
    public:
        [[nodiscard]] static ScalarType softmax(const VectorType& v, size_t label);
        [[nodiscard]] static ScalarType crossEntropy(const VectorType& v, size_t label);
        [[nodiscard]] static ScalarType polar_rate(ScalarType train_loss, ScalarType valid_loss);
        [[nodiscard]] static ScalarType mixed_loss(ScalarType train_loss, ScalarType valid_loss);
    };

    template<class ScalarType>
    ScalarType Loss<ScalarType>::softmax(const VectorType& v, size_t label) {
        const ScalarType maximum = v.max();
        const ScalarType factor = reciprocal(exp(v - maximum).sum());
        return exp(v.calc(label) - maximum) * factor;
    }

    template<class ScalarType>
    ScalarType Loss<ScalarType>::crossEntropy(const VectorType& v, size_t label) {
        assert(label < v.getLength() && "[Error]: The label is not exist");
        return -ln(softmax(v, label) + ScalarType(std::numeric_limits<ValueType>::min())); //Add minimum to avoid ln(0)
    }
    /**
     * \returns polarization rate, the lower the better, minus value means overfitting
     */
    template<class ScalarType>
    ScalarType Loss<ScalarType>::polar_rate(ScalarType train_loss, ScalarType valid_loss) {
        const ScalarType total = train_loss + valid_loss;
        const ScalarType delta = train_loss - valid_loss;
        if (abs(delta) * ScalarType(std::numeric_limits<ScalarType>::epsilon()) >= total)
            return ScalarType(0);
        return delta / total * ScalarType(2);
    }

    template<class ScalarType>
    ScalarType Loss<ScalarType>::mixed_loss(ScalarType train_loss, ScalarType valid_loss) {
        return std::max(train_loss, valid_loss) + abs(train_loss - valid_loss);
    }
}
