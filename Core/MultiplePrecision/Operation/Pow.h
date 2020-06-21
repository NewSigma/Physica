/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POW_H
#define PHYSICA_POW_H

#include "Core/Header/Scalar.h"
#include "Core/Header/SystemBits.h"
#include "Core/MultiplePrecision/Util/Bitwise.h"
//TODO Not tested
namespace Physica::Core {
    //Compute a ^ unit.
    inline Scalar powWord(const Scalar& a, ScalarUnit unit) {
        Scalar result(a);
        const auto lastUnitBits = countLeadingZeros(unit);
        for(int j = 0; j < ScalarUnitWidth - lastUnitBits; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    //Compute a ^ unit, the highest bit of unit must be set.
    inline Scalar powFullWord(const Scalar& a, ScalarUnit unit) {
        Scalar result(a);
        for(int j = 0; j < 64; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }

    inline Scalar powNumerical(const Scalar& a, const Scalar& n) {
        const auto size = n.getSize();
        Scalar result(a);
        if(n.getLength() < 0)
            result = reciprocal(a);

        for(int i = 0; i < size - 1; ++i)
            result = powFullWord(result, n[i]);
        result = powWord(result, n[size - 1]);
        return result;
    }
}

#endif