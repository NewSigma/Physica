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
#ifndef PHYSICA_SUBBASIC_H
#define PHYSICA_SUBBASIC_H

namespace Physica::Core {
    /*
     * Calculate n1 - n2, the result must be one word.
     */
    inline void sub2WordByWord(MPUnit& high, MPUnit& low
            , MPUnit n1_high, MPUnit n1_low, MPUnit n2) {
        low = n1_low - n2;
        high = n1_high - low > n1_low;
    }
    //Another version of sub2WordByWord()
    inline void sub2WordByWord(MPUnit& high, MPUnit& low, MPUnit n2) {
        //Use n2 as a temp variable, instead of creating a new temp variable.
        n2 = low - n2;
        high -= low > n2;
        low = n2;
    }
    /*
     * Calculate arr1 - arr2. len is the length of arr1 and arr2.
     * If arr1 >= arr2 the function will return 0, if not, the calculation is failed and return true.
     */
    inline bool subArrByArr(MPUnit* __restrict result, const MPUnit* __restrict arr1
            , const MPUnit* __restrict arr2, size_t len) {
        MPUnit carry = 0, temp1, temp2;
        bool pre_carry;
        for(size_t i = 0; i < len; ++i) {
            temp1 = arr1[i];
            temp2 = temp1 - arr2[i];
            pre_carry = temp1 < temp2;
            temp1 = temp2 - carry;
            result[i] = temp1;
            carry = pre_carry | (temp1 > temp2);
        }
        return carry;
    }
    //Another version of subArrByArr(), calculate arr1 -= arr2.
    inline bool subArrByArrEq(MPUnit* __restrict arr1, const MPUnit* __restrict arr2, size_t len) {
        MPUnit carry = 0, temp1, temp2;
        bool pre_carry;
        for(size_t i = 0; i < len; ++i) {
            temp1 = arr1[i];
            temp2 = temp1 - arr2[i];
            pre_carry = temp1 < temp2;
            temp1 = temp2 - carry;
            arr1[i] = temp1;
            carry = pre_carry | (temp1 > temp2);
        }
        return carry;
    }
}

#endif