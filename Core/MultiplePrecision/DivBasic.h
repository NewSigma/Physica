/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_DIVBASIC_H
#define PHYSICA_DIVBASIC_H

#include "Core/Header/SystemBits.h"
#include "MulBasic.h"
/*
 * Return the precomputed reciprocal.
 * Reference: T. Granlund and N. M¨oller, “Division of integers large and small”, to appear.
 */
inline NumericalUnit getInverse(NumericalUnit unit) {
#ifdef PHYSICA_64BIT
    unsigned long unit0 = unit & 1U;
    unsigned long unit9 = unit >> 55U;
    unsigned long unit40 = (unit >> 24U) + 1;
    unsigned long unit63 = (unit >> 1U) + unit0;
    unsigned long v0 = ((static_cast<unsigned long>(1) << 19U) - 3 * (static_cast<unsigned long>(1) << 8U)) / unit9;
    unsigned long v1 = (static_cast<unsigned long>(1) << 11U) * v0 - ((v0 * v0 * unit40) >> 40U) - 1;
    unsigned long v2 = (static_cast<unsigned long>(1) << 13U) * v1
                       + ((((static_cast<unsigned long>(1) << 60U) - v1 * unit40) * v1) >> 47U);
    unsigned long e = (v2 >> 1U) * unit0 - v2 * unit63;
    unsigned long v3 = (static_cast<unsigned long>(1) << 31U) * v2 + (mulWordByWordHigh(v2, e) >> 1U);
    unsigned long high;
    auto low = mulWordByWord(high, v3, unit);
    unsigned long v4 = v3 - unit - high - (low + unit < low);
    return v4;
#endif
#ifdef Physica_32BIT
#endif
}

#endif