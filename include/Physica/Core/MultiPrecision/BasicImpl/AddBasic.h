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
#ifndef PHYSICA_ADDBASIC_H
#define PHYSICA_ADDBASIC_H

namespace Physica::Core {
    /*!
     * \len is the length of \from. Length of \to should not less than \len. length of \result should be \len at least.
     */
    inline ScalarUnit addArrWithArr(ScalarUnit* __restrict result
            , const ScalarUnit* __restrict from, const ScalarUnit* __restrict to, size_t len) {
        ScalarUnit carry = 0, from_i, temp;
        for(int i = 0; i < len; ++i) {
            from_i = from[i];
            temp = from_i + to[i];
            result[i] = temp + carry;
            carry = temp < from_i;
        }
        return carry;
    }
    /*!
     * \len is the length of \from. Length of \to should not less than \len.
     */
    inline ScalarUnit addArrWithArrEq(const ScalarUnit* __restrict from, ScalarUnit* __restrict to, size_t len) {
        ScalarUnit carry = 0, from_i, temp;
        for(int i = 0; i < len; ++i) {
            from_i = from[i];
            temp = from_i + to[i];
            to[i] = temp + carry;
            carry = temp < from_i;
        }
        return carry;
    }
}

#endif