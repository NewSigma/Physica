/*
 * Copyright 2021-2022 WeiBo He.
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

#include "NumCharacter.h"

namespace Physica::Core {
    template<class VectorType>
    class LinearFit {
        using ScalarType = typename VectorType::ScalarType;
        using Param = std::pair<ScalarType, ScalarType>;
    public:
        static Param fit(const LValueVector<VectorType>& x, const LValueVector<VectorType>& y);
        static ScalarType relatedCoeff(const LValueVector<VectorType>& x, const LValueVector<VectorType>& y);
        static ScalarType deviation(const LValueVector<VectorType>& x, const LValueVector<VectorType>& y, Param param);
    };

    template<class VectorType>
    typename LinearFit<VectorType>::Param
    LinearFit<VectorType>::fit(const LValueVector<VectorType>& x, const LValueVector<VectorType>& y) {
        assert(x.getLength() == y.getLength());
        const size_t length = x.getLength();
        ScalarType xy_sum(0);
        ScalarType x_square_sum(0);
        ScalarType x_sum(0);
        ScalarType y_sum(0);
        for (size_t i = 0; i < length; ++i) {
            const auto& x_i = x[i];
            const auto& y_i = y[i];
            xy_sum += x_i * y_i;
            y_sum += y_i;
            x_square_sum += square(x_i);
            x_sum += x_i;
        }
        const ScalarType length_1 = reciprocal(ScalarType(length));
        const ScalarType numerator = xy_sum - x_sum * y_sum * length_1;
        const ScalarType denominator = x_square_sum - square(x_sum) * length_1;
        const ScalarType k = numerator / denominator;
        return {k, (y_sum - k * x_sum) * length_1};
    }

    template<class VectorType>
    typename LinearFit<VectorType>::ScalarType
    LinearFit<VectorType>::relatedCoeff(const LValueVector<VectorType>& x, const LValueVector<VectorType>& y) {
        return covariance(x, y) / sqrt(variance(x) * variance(y));
    }

    template<class VectorType>
    typename LinearFit<VectorType>::ScalarType
    LinearFit<VectorType>::deviation(const LValueVector<VectorType>& x, const LValueVector<VectorType>& y, Param param) {
        ScalarType result = square(y - param.first * x - param.second).sum();
        result = sqrt(result / ScalarType(x.getLength() - 2));
        return result;
    }
}
