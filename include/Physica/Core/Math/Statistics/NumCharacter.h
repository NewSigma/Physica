/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Utils/Container/Array/Array.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/LValueVector.h"

namespace Physica::Core {
    template<class VectorType>
    typename VectorType::ScalarType mean(const LValueVector<VectorType>& x) {
        using ScalarType = typename VectorType::ScalarType;
        return x.getDerived().sum() / ScalarType(x.getLength());
    }

    template<class VectorType>
    typename VectorType::ScalarType variance(const LValueVector<VectorType>& x) {
        using ScalarType = typename VectorType::ScalarType;
        const size_t length = x.getLength();
        const ScalarType x_mean = mean(x);
        return square(x - x_mean).sum() / ScalarType(length);
    }

    template<class VectorType>
    typename VectorType::ScalarType deviation(const LValueVector<VectorType>& x) {
        return sqrt(variance(x));
    }

    template<class VectorType>
    typename VectorType::ScalarType covariance(const LValueVector<VectorType>& x, const LValueVector<VectorType>& y) {
        assert(x.getLength() == y.getLength());
        using ScalarType = typename VectorType::ScalarType;
        const ScalarType x_mean = mean(x);
        const ScalarType y_mean = mean(y);
        return (x - x_mean) * (y - y_mean) / ScalarType(x.getLength());
    }

    template<class VectorType>
    typename VectorType::ScalarType skew(const LValueVector<VectorType>& x) {
        using ScalarType = typename VectorType::ScalarType;
        VectorType temp = x;
        const ScalarType factor = reciprocal(deviation(x));
        temp = (x - mean(x)) * factor;
        temp = hadamard(square(temp), temp);
        return mean(temp);
    }

    template<class VectorType>
    typename VectorType::ScalarType kurt(const LValueVector<VectorType>& x) {
        using ScalarType = typename VectorType::ScalarType;
        VectorType temp = x;
        const ScalarType factor = reciprocal(deviation(x));
        temp = (x - mean(x)) * factor;
        temp = square(square(temp));
        return mean(temp);
    }
}
