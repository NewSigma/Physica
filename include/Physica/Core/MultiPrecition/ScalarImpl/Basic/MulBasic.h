/*
 * Copyright 2019 WeiBo He.
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
#ifndef PHYSICA_MULBASIC_H
#define PHYSICA_MULBASIC_H

#include "Physica/SystemBits.h"

namespace Physica::Core {
    /*!
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
        unsigned long n1_low = n1 & ScalarUnitLowerMask;
        unsigned long n1_high = n1 >> (64U / 2U);
        unsigned long n2_low = n2 & ScalarUnitLowerMask;
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
    /*!
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
        ScalarUnit n1_low = n1 & ScalarUnitLowerMask;
        ScalarUnit n1_high = n1 >> (ScalarUnitWidth / 2U);
        ScalarUnit n2_low = n2 & ScalarUnitLowerMask;
        ScalarUnit n2_high = n2 >> (ScalarUnitWidth / 2U);

        auto ll = n1_low * n2_low;
        auto lh = n1_low * n2_high;
        auto hl = n1_high * n2_low;
        auto hh = n1_high * n2_high;

        lh += ll >> (ScalarUnitWidth / 2U);
        lh += hl;
        hh += static_cast<ScalarUnit>(lh < hl) << (ScalarUnitWidth / 2U);
        high = hh + (lh >> (ScalarUnitWidth / 2U));
        low = (lh << (ScalarUnitWidth / 2U)) + (ll & ScalarUnitLowerMask);
    #endif
    }
    /*!
     * Multiply the array @param arr with @param n. Write the result to array @param result.
     * Length of result should at least as long as arr.
     */
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
    /*!
     * Multiply the array @param arr with @param n. Add the result to array @param result.
     * Length of result should at least as long as arr.
     */
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
    /*!
     * Multiply the array @param arr with @param n. Subtract the result from array @param result.
     * Length of result should at least as long as arr.
     */
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