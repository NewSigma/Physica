/*
 * Copyright 2021 WeiBo He.
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
#include "Physica/Core/Math/Discrete/Combination.h"

namespace Physica::Core {
    /**
     * Calculate the number of arrangement $A_i1^i2$.
     * Using the definition, may be possible to optimize.
     */
    Integer arrangement(const Integer& i1, const Integer& i2) {
        assert(i1 > i2);
        const Integer critical = i1 - i2;
        Integer temp(i1);
        Integer result = 1;
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
    Integer combination(const Integer& i1, const Integer& i2) {
        assert(i1 > i2);
        const Integer critical = i1 - i2;
        Integer result = 1;
        const bool flag = critical > i2;
        const Integer& great = flag ? critical : i1;
        const Integer& small = flag ? i2 : critical;

        Integer temp(i1);
        while(temp > great) {
            result *= temp;
            --temp;
        }
        return result / factorial(small);
    }
}