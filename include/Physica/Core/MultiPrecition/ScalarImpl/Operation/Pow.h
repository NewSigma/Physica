/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POW_H
#define PHYSICA_POW_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/MultiPrecition/ScalarImpl/Util/Bitwise.h"
#include "Physica/Core/MultiPrecition/ElementaryFunction.h"
//TODO Not tested
namespace Physica::Core {
    //!Compute a ^ unit.
    inline Scalar<MultiPrecision, false> powWord(const Scalar<MultiPrecision, false>& a, ScalarUnit unit) {
        Scalar<MultiPrecision, false> result(a);
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
    inline Scalar<MultiPrecision, false> powFullWord(const Scalar<MultiPrecision, false>& a, ScalarUnit unit) {
        Scalar<MultiPrecision, false> result(a);
        for(int j = 0; j < 64; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }

    inline Scalar<MultiPrecision, false> powNumerical(
            const Scalar<MultiPrecision, false>& a, const Scalar<MultiPrecision, false>& n) {
        const auto size = n.getSize();
        Scalar<MultiPrecision, false> result(a);
        if(n.getLength() < 0)
            result = reciprocal(a);

        for(int i = 0; i < size - 1; ++i)
            result = powFullWord(result, n[i]);
        result = powWord(result, n[size - 1]);
        return result;
    }
}

#endif