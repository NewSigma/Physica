/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_BITWISE_H
#define PHYSICA_BITWISE_H

extern const NumericalUnit numericalUnitHighestBitMask;

//Possibly use asm to speed up.
inline unsigned int countLeadingZeros(NumericalUnit n) {
    if(n == 0)
        return NumericalUnitWidth;
    int count = 0;
    while((n & numericalUnitHighestBitMask) == 0) {
        ++count;
        n <<= 1U;
    }
    return count;
}

#endif