/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SUBBASIC_H
#define PHYSICA_SUBBASIC_H

#include "Core/Header/SystemBits.h"

namespace Physica::Core {
    /*
     * Calculate n1 - n2, the result must be one word.
     */
    inline void sub2WordByWord(NumericalUnit& high, NumericalUnit& low
            , NumericalUnit n1_high, NumericalUnit n1_low, NumericalUnit n2) {
        low = n1_low - n2;
        high = n1_high - low > n1_low;
    }
    //Another version of sub2WordByWord()
    inline void sub2WordByWord(NumericalUnit& high, NumericalUnit& low, NumericalUnit n2) {
        //Use n2 as a temp variable, instead of creating a new temp variable.
        n2 = low - n2;
        high -= low > n2;
        low = n2;
    }
    /*
     * Calculate arr1 - arr2. len is the length of arr1 and arr2.
     * If arr1 >= arr2 the function will return 0, if not, the calculation is failed and return true.
     */
    inline bool subArrByArr(NumericalUnit* result, const NumericalUnit* arr1, const NumericalUnit* arr2, size_t len) {
        NumericalUnit carry = 0, temp, pre_carry;
        for(int i = 0; i < len; ++i) {
            result[i] = arr1[i] - arr2[i];
            pre_carry = result[i] > arr1[i];
            temp = result[i] - carry;
            carry = pre_carry | (temp > result[i]);
            result[i] = temp;
        }
        return carry;
    }
    //Another version of subArrByArr()
    inline bool subArrByArrEq(NumericalUnit* arr1, const NumericalUnit* arr2, size_t len) {
        NumericalUnit carry = 0, temp, pre_carry;
        for(int i = 0; i < len; ++i) {
            temp = arr1[i] - arr2[i];
            pre_carry = temp > arr1[i];
            arr1[i] = temp;
            arr1[i] -= carry;
            carry = pre_carry | (arr1[i] > temp);
        }
        return carry;
    }
}

#endif