/*
 * Copyright 2020-2021 WeiBo He.
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

#include "IntegerArithmetic.h"

namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    Integer::Integer(const Scalar<option, errorTrack>& s) : Integer(static_cast<int>(s.getTrivial())) {}

    template<bool errorTrack>
    Integer::Integer(const Scalar<MultiPrecision, errorTrack>& s) {
        const auto power = s.getPower();
        if (power < 0) {
            byte = reinterpret_cast<MPUnit*>(malloc(sizeof(MPUnit)));
            byte[0] = 0;
            length = 1;
        }
        length = power + 1;
        const size_t size = length * sizeof(MPUnit);
        byte = reinterpret_cast<MPUnit*>(malloc(size));
        memcpy(byte, s.byte, length * sizeof(MPUnit));
    }
    /**
     * Returns true if i1 and i2 has the same sign. Both i1 and i2 do not equal to zero.
     * This function provide a quick sign check compare to using isPositive() and isNegative().
     */
    inline bool Integer::matchSign(const Integer& i1, const Integer& i2) {
        assert(!i1.isZero() && !i2.isZero());
        return (i1.length ^ i2.length) >= 0;
    }

    inline void swap(Integer& i1, Integer& i2) noexcept {
        i1.swap(i2);
    }
}

#include "ElementaryFunction.h"
