/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PROBABILITYFUNCTIONIMPL_H
#define PHYSICA_PROBABILITYFUNCTIONIMPL_H

#ifndef PHYSICA_PROBABILITYFUNCTION_H
    #include "../ProbabilityFunction.h"
#endif

namespace Physica::Core {
    //!Return a real number between 0 and 1.
    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> randomScalar() {
        return Scalar<maxPrecision, errorTrack>(double(random()) / RAND_MAX);
    }
    //!Return a real number lowerBound and upperBound.
    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> randomScalar(
            const Scalar<maxPrecision, errorTrack>& lowerBound,
            const Scalar<maxPrecision, errorTrack>& upperBound) {
        auto s = randomScalar<maxPrecision, errorTrack>();
        return s * (upperBound - lowerBound) + lowerBound;
    }

    template<size_t maxPrecision, bool errorTrack>
    Scalar<maxPrecision, errorTrack> floor(const Scalar<maxPrecision, errorTrack>& n) {
        if(n.isInteger())
            return Scalar(n);
        const auto size = n.getSize();
        const auto power = n.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = n.isNegative() ? -length : length;
        Scalar<maxPrecision, errorTrack> result(length, power);
        for(int i = 0; i < length; ++i)
            result[i] = n[i];
        return result;
    }

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> ceil(const Scalar<maxPrecision, errorTrack>& n) {
        return ++floor(n);
    }
}

#endif