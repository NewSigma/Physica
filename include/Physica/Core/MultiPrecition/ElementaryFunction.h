/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ELEMENTARYFUNCTION_H
#define PHYSICA_ELEMENTARYFUNCTION_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> square(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> reciprocal(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sqrt(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> factorial(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> ln(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    Scalar<type, errorTrack1 | errorTrack2> log(const Scalar<type, errorTrack1>& s, const Scalar<type, errorTrack2>& a);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> exp(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> cos(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sin(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> tan(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sec(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> csc(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> cot(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccos(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arcsin(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arctan(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arcsec(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccsc(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccot(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> cosh(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sinh(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> tanh(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sech(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> csch(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> coth(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccosh(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arcsinh(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arctanh(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arcsech(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccsch(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccoth(const Scalar<type, errorTrack>& s);
}

#include "FunctionsImpl/ElementaryFunctionImpl.h"

#endif