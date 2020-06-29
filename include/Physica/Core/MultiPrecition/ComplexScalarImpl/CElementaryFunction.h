/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CELEMENTARYFUNCTION_H
#define PHYSICA_CELEMENTARYFUNCTION_H
/*!
 * This file is part of implementations of \ComplexScalar.
 * Do not include this header file, include ComplexScalar.h instead.
 *
 * A number functions have not been implemented.
 */
namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> square(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> reciprocal(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sqrt(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> ln(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> exp(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> cos(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sin(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> tan(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sec(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> csc(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> cot(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> cosh(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sinh(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> tanh(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sech(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> csch(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> coth(const ComplexScalar<type, errorTrack>& c);
}

#include "FunctionImpl/CElementaryImpl.h"

#endif