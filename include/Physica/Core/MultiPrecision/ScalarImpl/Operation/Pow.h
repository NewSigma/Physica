/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_POW_H
#define PHYSICA_POW_H

namespace Physica::Core {
    //!Compute a ^ unit.
    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> powWord(
            const Scalar<MultiPrecision, errorTrack>& a, ScalarUnit unit) {
        Scalar<MultiPrecision, errorTrack> result(a);
        const auto lastUnitBits = countLeadingZeros(unit);
        for(unsigned int j = 0; j < ScalarUnitWidth - lastUnitBits; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    template<bool errorTrack>
    //!Compute a ^ unit, the highest bit of unit must be set.
    inline Scalar<MultiPrecision, errorTrack> powFullWord(
            const Scalar<MultiPrecision, errorTrack>& a, ScalarUnit unit) {
        Scalar<MultiPrecision, errorTrack> result(a);
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
     * Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009.45
     */
    template<bool errorTrack1, bool errorTrack2>
    inline Scalar<MultiPrecision, errorTrack1 | errorTrack2> powScalar(
            const Scalar<MultiPrecision, errorTrack1>& a, const Scalar<MultiPrecision, errorTrack2>& n) {
        const auto size = n.getSize();
        Scalar<MultiPrecision, errorTrack1 | errorTrack2> result(a);
        if(n.getLength() < 0)
            result = reciprocal(a);

        for(int i = 0; i < size - 1; ++i)
            result = powFullWord(result, n[i]);
        result = powWord(result, n[size - 1]);
        return result;
    }
}

#endif