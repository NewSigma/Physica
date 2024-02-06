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

#include "Loss.h"

namespace Physica::Core {
    template<class ScalarType>
    class Loss<device_obj<ScalarType>> {
        static_assert(ScalarType::isDifferentiable, "[Error]: Unexpected indifferentiable device scalar");
        using PlainScalar = typename ScalarType::PlainScalar;
        using DiffScalar = device_obj<Differentiable<PlainScalar, DiffMode::Reverse>>;
        using VectorType = device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>;
    public:
        [[nodiscard]] static DiffScalar crossEntropy(const VectorType& result, size_t label);
    };

    template<class ScalarType>
    typename Loss<device_obj<ScalarType>>::DiffScalar
    Loss<device_obj<ScalarType>>::crossEntropy(const VectorType& result, size_t label) {
        assert(label < result.getLength() && "[Error]: The label is not exist");
        return -ln(softmax(result).calc(label) + DiffScalar(std::numeric_limits<PlainScalar>::min())); //Add minimum to avoid ln(0)
    }
}
