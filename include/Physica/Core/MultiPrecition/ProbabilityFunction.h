/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PROBABILITYFUNCTION_H
#define PHYSICA_PROBABILITYFUNCTION_H

#include "Scalar.h"

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar();
    
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar(
            const Scalar<type, errorTrack>& lowerBound,
            const Scalar<type, errorTrack>& upperBound);

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> floor(const Scalar<type, errorTrack>& n);
    
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> ceil(const Scalar<type, errorTrack>& n);
}

#include "FunctionsImpl/ProbabilityFunctionImpl.h"

#endif