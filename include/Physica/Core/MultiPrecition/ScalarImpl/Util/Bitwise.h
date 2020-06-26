/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_BITWISE_H
#define PHYSICA_BITWISE_H

#include "Physica/Core/SystemBits.h"

namespace Physica::Core {
    //Possibly use asm to speed up.
    inline unsigned int countLeadingZeros(ScalarUnit n) {
        if(n == 0)
            return ScalarUnitWidth;
    #if UseASM
        ScalarUnit count;
        asm volatile (
                "bsrq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
        );
        (count) ^= 63U;
    #else
        ScalarUnit count = 0U;
        while((n & ScalarUnitHighestBitMask) == 0) {
            ++count;
            n <<= 1U;
        }
    #endif
        return count;
    }
}

#endif