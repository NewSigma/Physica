/*
 * Copyright 2024 Weibo He.
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
    /**
     * Returns true if s1 and s2 has the same sign. Both s1 and s2 do not equal to zero.
     * This function provide a quick sign check compare to using isPositive() and isNegative().
     */
    inline bool Scalar<MultiPrecision>::matchSign(const Scalar<MultiPrecision>& s1, const Scalar<MultiPrecision>& s2) {
        assert(!s1.isZero() && !s2.isZero());
        return (s1.length ^ s2.length) >= 0; //NOLINT Bitwise operator between two signed integer is intended.
    }
    /**
     * Cut zeros from the beginning.
     */
    inline void Scalar<MultiPrecision>::cutZero() {
        const int size = getSize();
        int id = size - 1;
        while(byte[id] == 0 && id > 0)
            --id;
        ++id;

        if(id != size) {
            int shorten = size - id;
            byte = reinterpret_cast<MPUnit*>(realloc(byte, id * sizeof(MPUnit)));
            length = length > 0 ? id : -id;
            auto temp = power;
            power = byte[id - 1] != 0 ? (temp - shorten) : 0;
        }
    }

    inline Scalar<MultiPrecision>& operator++(Scalar<MultiPrecision>& s) {
        s += BasicConst::getInstance()._1;
        return s;
    }

    inline Scalar<MultiPrecision>& operator--(Scalar<MultiPrecision>& s) {
        s -= BasicConst::getInstance()._1;
        return s;
    }

    inline Scalar<MultiPrecision> operator++(Scalar<MultiPrecision>& s, int) {
        Scalar<MultiPrecision> temp(s);
        s += BasicConst::getInstance()._1;
        return temp;
    }

    inline Scalar<MultiPrecision> operator--(Scalar<MultiPrecision>& s, int) {
        Scalar<MultiPrecision> temp(s);
        s -= BasicConst::getInstance()._1;
        return temp;
    }
}