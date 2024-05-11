/*
 * Copyright 2021-2024 WeiBo He.
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

#include <cassert>
#include <climits>

namespace Physica::Core {
    /**
     * Calculate the number of arrangement $A_i1^i2$.
     * Using the definition, may be possible to optimize.
     */
    template<class IntType>
    IntType arrangement(const IntType& i1, const IntType& i2) {
        static_assert(std::numeric_limits<IntType>::is_integer, "[Error]: Invalid IntType");
        assert(i1 > i2);
        const IntType critical = i1 - i2;
        IntType temp(i1);
        IntType result = 1;
        while(temp > critical) {
            result *= temp;
            --temp;
        }
        return result;
    }
    /**
     * Calculate the number of combination $C_i1^i2$.
     * Using the definition, may be possible to optimize.
     */
    template<class IntType>
    IntType combination(const IntType& i1, const IntType& i2) {
        static_assert(std::numeric_limits<IntType>::is_integer, "[Error]: Invalid IntType");
        assert(i1 > i2);
        const IntType critical = i1 - i2;
        IntType result = 1;
        const bool flag = critical > i2;
        const IntType& great = flag ? critical : i1;
        const IntType& small = flag ? i2 : critical;

        IntType temp(i1);
        while(temp > great) {
            result *= temp;
            --temp;
        }
        return result / factorial(small);
    }
}