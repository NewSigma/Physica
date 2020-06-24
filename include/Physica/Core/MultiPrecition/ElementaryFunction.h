/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ELEMENTARYFUNCTION_H
#define PHYSICA_ELEMENTARYFUNCTION_H

#include "Scalar.h"

namespace Physica::Core {
    template<size_t maxPrecision, bool errorTrack>
    Scalar<maxPrecision, errorTrack> square(const Scalar<maxPrecision, errorTrack>& n);

    template<size_t maxPrecision, bool errorTrack>
    inline Scalar<maxPrecision, errorTrack> reciprocal(const Scalar<maxPrecision, errorTrack>& n);

    Scalar sqrt_light(const Scalar& n);
    Scalar sqrt(const Scalar& n);
    Scalar factorial(const Scalar& n);
    Scalar ln_light(const Scalar& n);
    Scalar ln(const Scalar& n);
    Scalar log(const Scalar& n, const Scalar& a);
    Scalar exp(const Scalar& n);
    Scalar cos(const Scalar& n);
    Scalar sin(const Scalar& n);
    Scalar tan(const Scalar& n);
    Scalar sec(const Scalar& n);
    Scalar csc(const Scalar& n);
    Scalar cot(const Scalar& n);
    Scalar arccos(const Scalar& n);
    Scalar arcsin(const Scalar& n);
    Scalar arctan(const Scalar& n);
    Scalar arcsec(const Scalar& n);
    Scalar arccsc(const Scalar& n);
    Scalar arccot(const Scalar& n);
    Scalar cosh(const Scalar& n);
    Scalar sinh(const Scalar& n);
    Scalar tanh(const Scalar& n);
    Scalar sech(const Scalar& n);
    Scalar csch(const Scalar& n);
    Scalar coth(const Scalar& n);
    Scalar arccosh(const Scalar& n);
    Scalar arcsinh(const Scalar& n);
    Scalar arctanh(const Scalar& n);
    Scalar arcsech(const Scalar& n);
    Scalar arccsch(const Scalar& n);
    Scalar arccoth(const Scalar& n);
}

#include "FunctionsImpl/ElementaryFunctionImpl.h"

#endif