/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_BITWISE_H
#define PHYSICA_BITWISE_H

#include "Physica/SystemBits.h"

namespace Physica::Core {
    //Possibly use asm to speed up.
    inline unsigned int countLeadingZeros(ScalarUnit n) {
        if(n == 0)
            return ScalarUnitWidth;
        ScalarUnit count;
    #if UseASM
        asm volatile (
                "bsrq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
        );
        (count) ^= 63U;
    #else
        count = 0;
        while((n & ScalarUnitHighestBitMask) == 0) {
            ++count;
            n <<= 1U;
        }
    #endif
        return count;
    }

    inline unsigned int countBackZeros(unsigned long n) {
        if(n == 0)
            return ScalarUnitWidth;
        ScalarUnit count;
    #if UseASM
        asm volatile (
        "bsfq %1, %0\n\t"
        : "=r" (count)
        : "rm" (n)
        );
    #else
        count = 0;
        while((n & 1U) == 0) {
            ++count;
            n >>= 1U;
        }
    #endif
        return count;
    }
}

#endif