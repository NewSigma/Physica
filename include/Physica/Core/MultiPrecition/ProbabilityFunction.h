/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PROBABILITYFUNCTION_H
#define PHYSICA_PROBABILITYFUNCTION_H

#include "Scalar.h"

namespace Physica::Core {
    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> randomScalar();
    
    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> randomScalar(
            const Scalar<maxPrecision, errorTrack>& lowerBound,
            const Scalar<maxPrecision, errorTrack>& upperBound);

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> floor(const Scalar<maxPrecision, errorTrack>& n);
    
    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> ceil(const Scalar<maxPrecision, errorTrack>& n);
}

#include "FunctionsImpl/ProbabilityFunctionImpl.h"

#endif