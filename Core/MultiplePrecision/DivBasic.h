/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_DIVBASIC_H
#define PHYSICA_DIVBASIC_H

#include "Core/Header/SystemBits.h"
#include "MulBasic.h"
#include "qcompilerdetection.h"
/*
 * Return the precomputed reciprocal.
 * Reference: T. Granlund and N. M¨oller, “Division of integers large and small”, to appear.
 */
inline NumericalUnit getInverse(NumericalUnit unit) {
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
 * Reference: T. Granlund and N. M¨oller, “Division of integers large and small”, to appear.
 */
inline void div2WordByWord(NumericalUnit& quotient, NumericalUnit& remainder
        , NumericalUnit high, NumericalUnit low, NumericalUnit divisor) {
    NumericalUnit quotient2;
    mulWordByWord(quotient, quotient2, high, getInverse(divisor));

    auto temp = quotient2;
    quotient2 += low;
    quotient += high + (temp > quotient2) + 1;
    remainder = low - quotient * divisor;
    if(remainder > quotient2) {
        --quotient;
        remainder += divisor;
    }
    if(Q_UNLIKELY(remainder >= divisor)) {
        ++quotient;
        remainder -= divisor;
    }
}
/*
 * This is a simplified version of div2WordByWord(), which returns the quotient only.
 */
inline NumericalUnit div2WordByWordQ(NumericalUnit high, NumericalUnit low, NumericalUnit divisor) {
    NumericalUnit quotient, quotient2;
    mulWordByWord(quotient, quotient2, high, getInverse(divisor));
    quotient += high + 1;
    quotient2 += low;
    NumericalUnit remainder = low - quotient * divisor;
    if(remainder > quotient2) {
        --quotient;
        remainder += divisor;
    }
    if(Q_UNLIKELY(remainder >= divisor))
        ++quotient;
    return quotient;
}
/*
 * This is a simplified version of div2WordByWord(), which returns the remainder only.
 */
inline NumericalUnit div2WordByWordR(NumericalUnit high, NumericalUnit low, NumericalUnit divisor) {
    NumericalUnit quotient, quotient2;
    mulWordByWord(quotient, quotient2, high, getInverse(divisor));

    auto temp = quotient2;
    quotient2 += low;
    quotient += high + (temp > quotient2) + 1;
    auto remainder = low - quotient * divisor;
    if(remainder > quotient2)
        remainder += divisor;
    if(Q_UNLIKELY(remainder >= divisor))
        remainder -= divisor;
    return remainder;
}

#endif