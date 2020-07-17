/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMBERTHEORY_H
#define PHYSICA_NUMBERTHEORY_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    /*!
     * Return if s is a odd number.
     */
    template<ScalarType type>
    bool isOdd(const Scalar<type, false>& s) {
        return s.isInteger() && (static_cast<unsigned int>(s.getTrivial()) & 1U);
    }

    template<>
    bool isOdd(const Scalar<MultiPrecision, false>& s) {
        return s.isInteger() && (s[0] & 1U);
    }

    template<ScalarType type>
    bool isEven(const Scalar<type, false>& s) {
        return s.isInteger() && !(static_cast<unsigned int>(s.getTrivial()) & 1U);
    }

    template<>
    bool isEven(const Scalar<MultiPrecision, false>& s) {
        return s.isInteger() && !(s[0] & 1U);
    }
    /*!
     * Optimize: Make use of DP and the fact that $C_n^m = C_n^(n - m)$
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> bernoulli(const Scalar<type, errorTrack>& s) {
        Q_ASSERT(s.isInteger() && !s.isNegative());
        Scalar<type, errorTrack> result = Scalar<type, errorTrack>::getOne();
        if(s.isZero())
            return result;
        Scalar<type, false> temp = Scalar<type, false>::getZero();
        for(; temp < s; ++temp)
            result -= combination(s, temp) * bernoulli(temp) / (s - temp + 1);
        return result;
    }
}

#endif