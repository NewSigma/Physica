/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_MULBASIC_H
#define PHYSICA_MULBASIC_H

#include "Core/Header/SystemBits.h"

namespace Physica::Core {
    /*
     * This is simplified version of mulWordByWord(), which get the high Unit only.
     * It is slightly faster than mulWordByWord() if we are interested in the high Unit only.
     */
    //FIXME inline will cause incorrect results in release mode.
    inline ScalarUnit mulWordByWordHigh(ScalarUnit n1, ScalarUnit n2) {
    #if UseASM
        ScalarUnit result;
        #ifdef PHYSICA_64BIT
            asm (
                    "mulq %2"
                    : "=d"(result)
                    : "a"(n1), "rm"(n2)
            );
        #else
            asm (
                    "mull %2"
                    : "=d"(result)
                    : "a"(n1), "rm"(n2)
            );
        #endif
        return result;
    #else
        unsigned long n1_low = n1 & numericalUnitLowMask;
        unsigned long n1_high = n1 >> (64U / 2U);
        unsigned long n2_low = n2 & numericalUnitLowMask;
        unsigned long n2_high = n2 >> (64U / 2U);

        auto ll = n1_low * n2_low;
        auto lh = n1_low * n2_high;
        auto hl = n1_high * n2_low;
        auto hh = n1_high * n2_high;

        lh += ll >> (64U / 2U);
        lh += hl;
        hh += static_cast<unsigned long>(lh < hl) << (64U / 2U);

        return hh + (lh >> (64U / 2U));
    #endif
    }
    /*
     * On 64 bits machine(similar to 32 bit machine):
     * n1 * n2 = product(16 bytes) = carry(high 8 bytes) + ReturnValue(low bytes)
     */
    inline void mulWordByWord(ScalarUnit& high, ScalarUnit& low, ScalarUnit n1, ScalarUnit n2) {
    #if UseASM
        #ifdef PHYSICA_64BIT
            asm (
                    "mulq %3"
                    : "=d"(high), "=a"(low)
                    : "a"(n1), "rm"(n2)
            );
        #else
            asm (
                    "mull %3"
                    : "=d"(high), "=a"(low)
                    : "a"(n1), "rm"(n2)
            );
        #endif
    #else
        ScalarUnit n1_low = n1 & numericalUnitLowMask;
        ScalarUnit n1_high = n1 >> (ScalarUnitWidth / 2U);
        ScalarUnit n2_low = n2 & numericalUnitLowMask;
        ScalarUnit n2_high = n2 >> (ScalarUnitWidth / 2U);

        auto ll = n1_low * n2_low;
        auto lh = n1_low * n2_high;
        auto hl = n1_high * n2_low;
        auto hh = n1_high * n2_high;

        lh += ll >> (ScalarUnitWidth / 2U);
        lh += hl;
        hh += static_cast<ScalarUnit>(lh < hl) << (ScalarUnitWidth / 2U);
        high = hh + (lh >> (ScalarUnitWidth / 2U));
        low = (lh << (ScalarUnitWidth / 2U)) + (ll & numericalUnitLowMask);
    #endif
    }
    //Length of result should at least as long as arr.
    inline ScalarUnit mulArrByWord(ScalarUnit* __restrict result, const ScalarUnit* __restrict arr
            , size_t length, ScalarUnit n) {
        ScalarUnit carry = 0, high, low;
        for(size_t i = 0; i < length; ++i) {
            mulWordByWord(high, low, arr[i], n);
            low += carry;
            carry = (low < carry) + high;
            result[i] = low;
        }
        return carry;
    }
    //Length of result should at least as long as arr.
    inline ScalarUnit mulAddArrByWord(ScalarUnit* __restrict result, const ScalarUnit* __restrict arr
            , size_t length, ScalarUnit n) {
        ScalarUnit carry = 0, high, low;
        for(size_t i = 0; i < length; ++i) {
            mulWordByWord(high, low, arr[i], n);
            low += carry;
            carry = (low < carry) + high;
            result[i] += low;
            carry += result[i] < low;
        }
        return carry;
    }
    //Length of result should at least as long as arr.
    inline ScalarUnit mulSubArrByWord(ScalarUnit* __restrict result, const ScalarUnit* __restrict arr, size_t length, ScalarUnit n) {
        ScalarUnit carry = 0, high, low;
        for(size_t i = 0; i < length; ++i) {
            mulWordByWord(high, low, arr[i], n);
            low = result[i] - low;
            high += result[i] < low;
            result[i] = low;
            low -= carry;
            carry = high + (result[i] < low);
            result[i] = low;
        }
        return carry;
    }
}

#endif