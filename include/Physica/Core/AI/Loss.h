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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    template<class ScalarType>
    class Loss {
        using PlainScalar = typename ScalarType::PlainScalar;
        using VectorType = Vector<ScalarType>;
    public:
        [[nodiscard]] static ScalarType crossEntropy(const VectorType& result, const VectorType& answer);
        [[nodiscard]] static ScalarType polar_rate(ScalarType train_loss, ScalarType valid_loss);
        [[nodiscard]] static ScalarType mixed_loss(ScalarType train_loss, ScalarType valid_loss);
    };

    template<class ScalarType>
    ScalarType Loss<ScalarType>::crossEntropy(const VectorType& result, const VectorType& answer) {
        for (size_t i = 0; i < answer.getLength(); ++i)
            if (!answer[i].isZero())
                return -ln(softmax(result).calc(i) + ScalarType(std::numeric_limits<PlainScalar>::min())); //Add minimum to avoid ln(0)
        assert(false && "[Error]: Invalid argument");
        return ScalarType(0);
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
