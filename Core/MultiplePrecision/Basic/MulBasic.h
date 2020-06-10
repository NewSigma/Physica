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
    inline NumericalUnit mulWordByWordHigh(NumericalUnit n1, NumericalUnit n2) {
    #if UseASM
        NumericalUnit result;
        #ifdef PHYSICA_64BIT
            asm (
                    "movq %1, %%rax\n\t"
                    "mulq %2"
                    : "=d"(result)
                    : "rm"(n1), "r"(n2)
                    : "%rax"
            );
        #else
            asm (
                    "movl %1, %%rax\n\t"
                    "mull %2"
                    : "=d"(result)
                    : "rm"(n1), "rm"(n2)
                    : "%rax"
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
    inline void mulWordByWord(NumericalUnit& high, NumericalUnit& low, NumericalUnit n1, NumericalUnit n2) {
    #if UseASM
        #ifdef PHYSICA_64BIT
            asm (
                    "movq %2, %%rax\n\t"
                    "mulq %3"
                    : "=d"(high), "=a"(low)
                    : "rm"(n1), "r"(n2)
            );
        #else
            asm (
                    "movl %2, %%rax\n\t"
                    "mull %3"
                    : "=d"(high), "=a"(low)
                    : "rm"(n1), "rm"(n2)
            );
        #endif
    #else
        NumericalUnit n1_low = n1 & numericalUnitLowMask;
        NumericalUnit n1_high = n1 >> (NumericalUnitWidth / 2U);
        NumericalUnit n2_low = n2 & numericalUnitLowMask;
        NumericalUnit n2_high = n2 >> (NumericalUnitWidth / 2U);

        auto ll = n1_low * n2_low;
        auto lh = n1_low * n2_high;
        auto hl = n1_high * n2_low;
        auto hh = n1_high * n2_high;

        lh += ll >> (NumericalUnitWidth / 2U);
        lh += hl;
        hh += static_cast<NumericalUnit>(lh < hl) << (NumericalUnitWidth / 2U);
        high = hh + (lh >> (NumericalUnitWidth / 2U));
        low = (lh << (NumericalUnitWidth / 2U)) + (ll & numericalUnitLowMask);
    #endif
    }
    //Length of result should at least as long as arr.
    inline NumericalUnit mulArrByWord(NumericalUnit* result, const NumericalUnit* arr, size_t length, NumericalUnit n) {
        NumericalUnit carry = 0, high, low;
        for(size_t i = 0; i < length; ++i) {
            mulWordByWord(high, low, arr[i], n);
            low += carry;
            carry = (low < carry) + high;
            result[i] = low;
        }
        return carry;
    }
    //Length of result should at least as long as arr.
    inline NumericalUnit mulAddArrByWord(NumericalUnit* result, const NumericalUnit* arr, size_t length, NumericalUnit n) {
        NumericalUnit carry = 0, high, low;
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
    inline NumericalUnit mulSubArrByWord(NumericalUnit* result, const NumericalUnit* arr, size_t length, NumericalUnit n) {
        NumericalUnit carry = 0, high, low;
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