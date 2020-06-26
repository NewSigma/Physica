/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PROBABILITYFUNCTIONIMPL_H
#define PHYSICA_PROBABILITYFUNCTIONIMPL_H
/*!
 * This file is part of implementations of \ProbabilityFunction.
 * Do not include this header file, include ProbabilityFunction.h instead.
 */
namespace Physica::Core {
    //!Return a real number between 0 and 1.
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar() {
        return Scalar<type, errorTrack>(static_cast<double>(random()) / RAND_MAX);
    }
    //!Return a real number lowerBound and upperBound.
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar(
            const Scalar<type, errorTrack>& lowerBound,
            const Scalar<type, errorTrack>& upperBound) {
        auto castUpper = static_cast<Scalar<type, false>>(upperBound);
        auto castLower = static_cast<Scalar<type, false>>(lowerBound);
        auto s = randomScalar<type, errorTrack>();
        return s * (castUpper - castLower) + castLower;
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> floor(const Scalar<MultiPrecision, errorTrack>& n) {
        if(n.isInteger())
            return Scalar(n);
        const auto size = n.getSize();
        const auto power = n.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = n.isNegative() ? -length : length;
        Scalar<MultiPrecision, errorTrack> result(length, power);
        for(int i = 0; i < length; ++i)
            result[i] = n[i];
        return result;
    }
    /*!
     * Specialization for float and double.
     * Fix: May overflow.
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> floor(const Scalar<type, errorTrack>& n) {
        return Scalar<type, errorTrack>(static_cast<size_t>(n.getTrivial()));
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> ceil(const Scalar<type, errorTrack>& n) {
        return ++floor(n);
    }
}

#endif