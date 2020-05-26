/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_BITWISE_H
#define PHYSICA_BITWISE_H

namespace Physica::Core {
    //Possibly use asm to speed up.
    inline unsigned int countLeadingZeros(NumericalUnit n) {
        if(n == 0)
            return NumericalUnitWidth;
    #if UseASM
        NumericalUnit count;
        asm volatile (
                "bsrq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
        );
        (count) ^= 63U;
    #else
        NumericalUnit count = 0U;
        while((n & numericalUnitHighestBitMask) == 0) {
            ++count;
            n <<= 1U;
        }
    #endif
        return count;
    }
}

#endif