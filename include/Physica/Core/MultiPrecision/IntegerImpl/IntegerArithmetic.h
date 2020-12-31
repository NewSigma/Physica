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
#ifndef PHYSICA_INTEGERARITHMETIC_H
#define PHYSICA_SCALARARITHMETIC_H

#include <cassert>
#include "Physica/Core/MultiPrecision/BasicImpl/AddBasic.h"
#include "Physica/Core/MultiPrecision/BasicImpl/DivBasic.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Util/ArraySupport.h"
#include "Physica/Core/MultiPrecision/BasicImpl/Util/Bitwise.h"
#include "Physica/Core/Exception/DivideByZeroException.h"

namespace Physica::Core {
    /**
     * Both i1 and i2 must be positive integers.
     */
    inline Integer Integer::integerAddImpl(const Integer& i1, const Integer& i2) {
        assert(i1.isPositive() && i2.isPositive());
        const Integer* __restrict largeInt = i1.length > i2.length ? &i1 : &i2;
        const Integer* __restrict smallInt = i1.length > i2.length ? &i2 : &i1;
        int length = largeInt->length;
        auto* __restrict byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));
        memcpy(byte, largeInt->byte, length * sizeof(MPUnit));

        bool carry = addArrWithArrEq(smallInt->byte, byte, smallInt->length);
        //usableSmall is the index that we will add carry to.
        int carryToIndex = smallInt->length;
        while (carry != 0 && carryToIndex < length) {
            MPUnit temp = byte[carryToIndex] + 1;
            byte[carryToIndex] = temp;
            carry = temp < carry;
            ++carryToIndex;
        }
        if (carry) {
            ++length;
            byte = reinterpret_cast<MPUnit*>(realloc(byte, length * sizeof(MPUnit)));
            byte[length - 1] = 1;
        }
        return Integer(byte, length);
    }
    /**
     * Both i1 and i2 must be positive integers.
     */
    inline Integer Integer::integerSubImpl(const Integer& i1, const Integer& i2) {
        assert(i1.isPositive() && i2.isPositive());
        bool changeSign = i1.length <= i2.length;
        const Integer* __restrict largeInt = i1.length > i2.length ? &i1 : &i2;
        const Integer* __restrict smallInt = i1.length > i2.length ? &i2 : &i1;
        if (i1.length == i2.length && !absCompare(i1, i2)) {
            changeSign = true;
            std::swap(largeInt, smallInt);
        }
        int length = largeInt->length;
        auto* __restrict byte = reinterpret_cast<MPUnit*>(malloc(length * sizeof(MPUnit)));
        memcpy(byte, largeInt->byte, length * sizeof(MPUnit));

        bool carry = subArrByArrEq(byte, smallInt->byte, smallInt->length);
        int carryToIndex = smallInt->length;
        while(carry != 0) {
            assert(carryToIndex < length);
            MPUnit temp1 = byte[carryToIndex];
            MPUnit temp2 = temp1 - 1;
            byte[carryToIndex] = temp2;
            carry = temp1 < temp2;
            ++carryToIndex;
        }
        Integer result(byte, changeSign ? -length : length);
        result.cutZero();
        return result;
    }
}

#endif
