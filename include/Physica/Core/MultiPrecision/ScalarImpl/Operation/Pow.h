/*
 * Copyright 2020-2023 WeiBo He.
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

namespace Physica::Core {
    //!Compute a ^ unit.
    inline Scalar<MultiPrecision> powWord(const Scalar<MultiPrecision>& a, MPUnit unit) {
        Scalar<MultiPrecision> result(a);
        const auto lastUnitBits = countLeadingZeros(unit);
        for(unsigned int j = 0; j < MPUnitWidth - lastUnitBits; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }

    //!Compute a ^ unit, the highest bit of unit must be set.
    inline Scalar<MultiPrecision> powFullWord(const Scalar<MultiPrecision>& a, MPUnit unit) {
        Scalar<MultiPrecision> result(a);
        for(int j = 0; j < 64; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    /*!
     * Calculate a^n.
     *
     * Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009:45
     */
    inline Scalar<MultiPrecision> powScalar(const Scalar<MultiPrecision>& a, const Scalar<MultiPrecision>& n) {
        const auto size = n.getSize();
        Scalar<MultiPrecision> result(a);
        if(n.getLength() < 0)
            result = reciprocal(a);

        for(int i = 0; i < size - 1; ++i)
            result = powFullWord(result, n[i]);
        result = powWord(result, n[size - 1]);
        return result;
    }
}
