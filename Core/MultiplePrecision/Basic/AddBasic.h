/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ADDBASIC_H
#define PHYSICA_ADDBASIC_H

#include <cstddef>
#include "Core/Header/SystemBits.h"

namespace Physica::Core {
    /*
     * Both of the arrays have the same length. Allocated length of result should not less than len + 1.
     */
    inline void addArrWithArr(NumericalUnit* result, const NumericalUnit* from, const NumericalUnit* to, size_t len) {
        result[0] = 0;
        for(int i = 0; i < len; ++i) {
            NumericalUnit temp = from[i] + to[i];
            result[i] += temp;
            result[i + 1] = temp < from[i];
        }
    }
    /*
     * Both of the arrays have the same length. Allocated length of to should not less than len + 1.
     */
    inline NumericalUnit addArrWithArrEq(const NumericalUnit* from, NumericalUnit* to, size_t len) {
        NumericalUnit carry = 0;
        for(int i = 0; i < len; ++i) {
            NumericalUnit temp = from[i] + to[i];
            to[i] += temp + carry;
            carry = temp < from[i];
        }
        return carry;
    }
}

#endif