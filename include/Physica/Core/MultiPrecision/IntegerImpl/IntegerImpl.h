/*
 * Copyright 2020 WeiBo He.
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
#include "IntegerArithmetic.h"

namespace Physica::Core {
    /**
     * Returns true if i1 and i2 has the same sign. Both i1 and i2 do not equal to zero.
     * This function provide a quick sign check compare to using isPositive() and isNegative().
     */
    bool Integer::matchSign(const Integer& i1, const Integer& i2) {
        assert(!i1.isZero() && !i2.isZero());
        return (i1.length ^ i2.length) >= 0;
    }
}