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
    inline void sub2WordByWord(ScalarUnit& high, ScalarUnit& low
            , ScalarUnit n1_high, ScalarUnit n1_low, ScalarUnit n2) {
        low = n1_low - n2;
        high = n1_high - low > n1_low;
    }
    //Another version of sub2WordByWord()
    inline void sub2WordByWord(ScalarUnit& high, ScalarUnit& low, ScalarUnit n2) {
        //Use n2 as a temp variable, instead of creating a new temp variable.
        n2 = low - n2;
        high -= low > n2;
        low = n2;
    }
    /*
     * Calculate arr1 - arr2. len is the length of arr1 and arr2.
     * If arr1 >= arr2 the function will return 0, if not, the calculation is failed and return true.
     */
    inline bool subArrByArr(ScalarUnit* __restrict result, const ScalarUnit* __restrict arr1
            , const ScalarUnit* __restrict arr2, size_t len) {
        ScalarUnit carry = 0, temp1, temp2;
        bool pre_carry;
        for(int i = 0; i < len; ++i) {
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
    inline bool subArrByArrEq(ScalarUnit* __restrict arr1, const ScalarUnit* __restrict arr2, size_t len) {
        ScalarUnit carry = 0, temp1, temp2;
        bool pre_carry;
        for(int i = 0; i < len; ++i) {
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