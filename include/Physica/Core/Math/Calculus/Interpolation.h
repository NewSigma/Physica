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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/LValueVector.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009:81-82
     */
    template<class VectorType>
    typename VectorType::ScalarType lagrange(
            const LValueVector<VectorType>& x0,
            const LValueVector<VectorType>& y0,
            typename VectorType::ScalarType x) {
        using ScalarType = typename VectorType::ScalarType;
        assert(x0.getLength() == y0.getLength());

        const size_t length = x0.getLength();
        VectorType c = y0;
        VectorType d = y0;
        size_t ns = 0;
        /* Init ns */ {
            ScalarType delta = abs(x - x0[0]);
            for (size_t i = 1; i < length; ++i) {
                const ScalarType temp = abs(x - x0[i]);
                if (temp < delta) {
                    ns = i;
                    delta = temp;
                }
            }
        }
        ScalarType result = y0[ns--];
        for (size_t i = 1; i < length; ++i) {
            const size_t max_j = length - i;
            for (size_t j = 0; j < max_j; ++j) {
                const ScalarType delta1 = x0[j] - x;
                const ScalarType delta2 = x0[j + i] - x;
                const ScalarType factor = (c[j + 1] - d[j]) / (delta1 - delta2);
                c[j] = delta1 * factor;
                d[j] = delta2 * factor;
            }
            result += (2 * (ns + 1) < max_j) ? c[ns + 1] : d[ns--];
        }
        return result;
    }
}
