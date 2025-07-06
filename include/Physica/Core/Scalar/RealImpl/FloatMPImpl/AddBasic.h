/*
 * Copyright 2020-2025 Weibo He.
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

#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    /*!
     * \len is the length of \from. Length of \to should not less than \len. length of \result should be \len at least.
     */
    inline MPUnit addArrWithArr(MPUnit* __restrict result, const MPUnit* __restrict from, const MPUnit* __restrict to, size_t len) noexcept {
        MPUnit carry = 0;
        for (size_t i = 0; i < len; ++i) {
            MPUnit from_i = from[i];
            MPUnit temp = from_i + to[i] + carry;
            result[i] = temp;
            carry = temp < from_i;
        }
        return carry;
    }
    /*!
     * \len is the length of \from. Length of \to should not less than \len.
     */
    inline MPUnit addArrWithArrEq(const MPUnit* __restrict from, MPUnit* __restrict to, size_t len) noexcept {
        MPUnit carry = 0;
        for (size_t i = 0; i < len; ++i) {
            MPUnit from_i = from[i];
            MPUnit add = from_i + to[i] + carry;
            to[i] = add;
            carry = add < from_i;
        }
        return carry;
    }
}
