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
        return square(x - x_mean).sum() / ScalarType(length - 1);
    }

    template<class VectorType>
    typename VectorType::ScalarType variance(const LValueVector<VectorType>& x, typename VectorType::ScalarType prior_mean) {
        using ScalarType = typename VectorType::ScalarType;
        const size_t length = x.getLength();
        return square(x - prior_mean).sum() / ScalarType(length);
    }

    template<class VectorType>
    typename VectorType::ScalarType deviation(const LValueVector<VectorType>& x) {
        return sqrt(variance(x));
    }

    template<class VectorType>
    VectorType normalize(const LValueVector<VectorType>& x) {
        using ScalarType = typename VectorType::ScalarType;
        const ScalarType x_mean = mean(x);
        const ScalarType factor = reciprocal(deviation(x));
        return VectorType((x - x_mean) * factor);
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
        VectorType temp = normalize(x);
        temp = hadamard(square(temp), temp);
        const size_t length = x.getLength();
        const ScalarType factor = ScalarType(length * length) / ScalarType((length - 1) * (length - 2));
        return mean(temp) * factor;
    }

    template<class VectorType>
    typename VectorType::ScalarType kurt(const LValueVector<VectorType>& x) {
        using ScalarType = typename VectorType::ScalarType;
        VectorType temp = normalize(x);
        temp = square(temp);
        const ScalarType mean2 = mean(temp);
        temp = square(temp);
        const ScalarType mean1 = mean(temp);

        const size_t length = x.getLength();
        const ScalarType factor2 = ScalarType(length * length * 3) / ScalarType((length - 2) * (length - 3));
        const ScalarType factor1 = ScalarType(length * length * (length + 1)) / ScalarType((length - 1) * (length - 2) * (length - 3));
        return factor1 * mean1 - factor2 * mean2;
    }
}
