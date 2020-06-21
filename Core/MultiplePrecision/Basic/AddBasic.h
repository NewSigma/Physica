/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ADDBASIC_H
#define PHYSICA_ADDBASIC_H

#include <cstddef>
#include "Core/Header/SystemBits.h"

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