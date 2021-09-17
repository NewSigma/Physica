/*
 * Copyright 2019-2021 WeiBo He.
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

#include "SubBasic.h"
#include "MulBasic.h"

namespace Physica::Core {
    /*
     * Return the precomputed reciprocal.
     * Reference: T. Granlund and N. M¨oller, “Division of integers large and small”, to appear.
     */
    inline MPUnit getInverse(MPUnit unit) {
        unsigned long unit0 = unit & 1U;
    #ifdef PHYSICA_64BIT
        unsigned long unit9 = unit >> 55U;
        unsigned long unit40 = (unit >> 24U) + 1;
        unsigned long unit63 = (unit >> 1U) + unit0;
        unsigned long v0 = ((static_cast<unsigned long>(1) << 19U) - 3 * (static_cast<unsigned long>(1) << 8U)) / unit9;
        unsigned long v1 = (static_cast<unsigned long>(1) << 11U) * v0 - ((v0 * v0 * unit40) >> 40U) - 1;
        unsigned long v2 = (static_cast<unsigned long>(1) << 13U) * v1
                           + ((((static_cast<unsigned long>(1) << 60U) - v1 * unit40) * v1) >> 47U);
        unsigned long e = (v2 >> 1U) * unit0 - v2 * unit63;
        unsigned long v3 = (static_cast<unsigned long>(1) << 31U) * v2 + (mulWordByWordHigh(v2, e) >> 1U);
        unsigned long high, low;
        mulWordByWord(high, low, v3, unit);
        unsigned long v4 = v3 - unit - high - (low + unit < low);
        return v4;
    #endif
    #ifdef PHYSICA_32BIT
        unsigned long unit10 = unit >> 22U;
        unsigned long unit21 = (unit >> 11U) + 1;
        unsigned long unit31 = (unit >> 1U) + unit10;
        unsigned long v0 = ((static_cast<unsigned long>(1) << 24U) - (static_cast<unsigned long>(1) << 14U)
                + (static_cast<unsigned long>(1) << 9U)) / unit10;
        unsigned long v1 = (static_cast<unsigned long>(1) << 4U) * v0 - mulWordByWordHigh(v0 * v0, unit21) - 1;
        unsigned long e = (v1 >> 1U) * unit0 - v1 * unit31;
        unsigned long v2 = (static_cast<unsigned long>(1) << 15U) * v1 - (mulWordByWordHigh(v1, e) >> 1U);
        unsigned long high, low;
        mulWordByWord(high, low, v2, unit);
        unsigned long v3 = v2 - unit - high - (low + unit < low);
        return v3;
    #endif
    }
    /*
     * Calculate (high, low) / divisor.
     * Assume high < divisor and divisor >= 2^(__WORDSIZE - 1).
     *
     * A full word here indicates that the highest bit of divisor is set.
     *
     * Reference: T. Granlund and N. M¨oller, “Division of integers large and small”, to appear.
     */
    inline void div2WordByFullWord(MPUnit& quotient, MPUnit& remainder
            , MPUnit high, MPUnit low, MPUnit divisor) {
        assert(high < divisor);
        assert(divisor & (static_cast<MPUnit>(1) << (PhysicaWordSize - 1U)));
        MPUnit quotient2;
        mulWordByWord(quotient, quotient2, high, getInverse(divisor));

        auto temp = quotient2;
        quotient2 += low;
        quotient += high + (temp > quotient2) + 1;
        remainder = low - quotient * divisor;
        if(remainder > quotient2) {
            --quotient;
            remainder += divisor;
        }
        if(remainder >= divisor) {
            ++quotient;
            remainder -= divisor;
        }
    }
    /*
     * This is a simplified version of div2WordByFullWord(), which returns the quotient only.
     */
    inline MPUnit div2WordByFullWordQ(MPUnit high, MPUnit low, MPUnit divisor) {
        assert(high < divisor);
        assert(divisor & (static_cast<MPUnit>(1) << (PhysicaWordSize - 1U)));
        MPUnit quotient, quotient2;
        mulWordByWord(quotient, quotient2, high, getInverse(divisor));

        auto temp = quotient2;
        quotient2 += low;
        quotient += high + (temp > quotient2) + 1;
        MPUnit remainder = low - quotient * divisor;
        if(remainder > quotient2) {
            --quotient;
            remainder += divisor;
        }
        if(remainder >= divisor)
            ++quotient;
        return quotient;
    }
    /*
     * This is a simplified version of div2WordByFullWord(), which returns the remainder only.
     */
    inline MPUnit div2WordByFullWordR(MPUnit high, MPUnit low, MPUnit divisor) {
        assert(high < divisor);
        assert(divisor & (static_cast<MPUnit>(1) << (PhysicaWordSize - 1U)));
        MPUnit quotient, quotient2;
        mulWordByWord(quotient, quotient2, high, getInverse(divisor));

        auto temp = quotient2;
        quotient2 += low;
        quotient += high + (temp > quotient2) + 1;
        auto remainder = low - quotient * divisor;
        if(remainder > quotient2)
            remainder += divisor;
        if(remainder >= divisor)
            remainder -= divisor;
        return remainder;
    }
    /*
     * Acquirement:
     * The quotient must be one word.
     * len is the length of divisor, length of dividend should equals to len + 1.
     * FullArr here indicates that the highest bit of divisor is set.
     * dividend[len] < divisor[len - 1].
     * len >= 1 to avoid invalid read.
     *
     * Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009:4-8
     */
    inline MPUnit divArrByFullArrWith1Word(const MPUnit* __restrict dividend
            , const MPUnit* __restrict divisor, size_t len) {
        assert(len >= 1);
        MPUnit q = dividend[len] >= divisor[len - 1] ? MPUnitMax :
                       div2WordByFullWordQ(dividend[len], dividend[len - 1], divisor[len - 1]);
        if(len == 1) //May be ask len > 1 to avoid branches.
            return q;
        MPUnit temp[2]{dividend[len - 1], dividend[len]};
        /* Calculate temp - q * divisor[len - 1] */ {
            MPUnit temp_1[2];
            mulWordByWord(temp_1[1], temp_1[0], q, divisor[len - 1]);
            subArrByArrEq(temp, temp_1, 2);
        }
        MPUnit q_divisor_high, q_divisor_low;
        mulWordByWord(q_divisor_high, q_divisor_low, q, divisor[len - 2]);
        sub2WordByWord(q_divisor_high, q_divisor_low, dividend[len - 2]);
        if(temp[1] == 0 && (temp[0] < q_divisor_high || (temp[0] == q_divisor_high && q_divisor_low == 0)))
            --q;

        MPUnit* __restrict n = new MPUnit[len + 1];
        n[len] = mulArrByWord(n, divisor, len, q);
        //Judge dividend < n or not. If dividend < n, we have to do carry.
        bool carry = true;
        for(size_t i = len; i <= len; --i) { //i <= len makes use of overflow
            if (dividend[i] > n[i]) {
                carry = false;
                break;
            }
        }
        delete[] n;
        return carry ? (q - 1) : q;
    }
}
