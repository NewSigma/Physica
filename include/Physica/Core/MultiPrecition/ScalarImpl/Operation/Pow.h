/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POW_H
#define PHYSICA_POW_H

//TODO Not tested
namespace Physica::Core {
    //!Compute a ^ unit.
    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> powWord(
            const Scalar<MultiPrecision, errorTrack>& a, ScalarUnit unit) {
        Scalar<MultiPrecision, errorTrack> result(a);
        const auto lastUnitBits = countLeadingZeros(unit);
        for(int j = 0; j < ScalarUnitWidth - lastUnitBits; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    template<bool errorTrack>
    //!Compute a ^ unit, the highest bit of unit must be set.
    inline Scalar<MultiPrecision, errorTrack> powFullWord(
            const Scalar<MultiPrecision, errorTrack>& a, ScalarUnit unit) {
        Scalar<MultiPrecision, errorTrack> result(a);
        for(int j = 0; j < 64; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    /*!
     * Calculate a^n.
     *
     * Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009.45
     */
    template<bool errorTrack1, bool errorTrack2>
    inline Scalar<MultiPrecision, errorTrack1 | errorTrack2> powScalar(
            const Scalar<MultiPrecision, errorTrack1>& a, const Scalar<MultiPrecision, errorTrack2>& n) {
        const auto size = n.getSize();
        Scalar<MultiPrecision, errorTrack1 | errorTrack2> result(a);
        if(n.getLength() < 0)
            result = reciprocal(a);

        for(int i = 0; i < size - 1; ++i)
            result = powFullWord(result, n[i]);
        result = powWord(result, n[size - 1]);
        return result;
    }
}

#endif