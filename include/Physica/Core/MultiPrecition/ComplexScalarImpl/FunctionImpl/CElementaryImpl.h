/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CELEMENTARYIMPL_H
#define PHYSICA_CELEMENTARYIMPL_H
/*!
 * This file is part of implementations of \ComplexScalar.
 * Do not include this header file, include ComplexScalar.h instead.
 */
namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> square(const ComplexScalar<type, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<type, errorTrack>(square(real) - square(imagine), (real * imagine) << 1);
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> reciprocal(const ComplexScalar<type, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        const auto divisor = reciprocal(square(real) + square(imagine));
        return ComplexScalar<type, errorTrack>(real * divisor, -imagine * divisor);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sqrt(const ComplexScalar<type, errorTrack>& c) {
        const auto n = norm(c);
        const auto a = arg(c) >> 1;
        return ComplexScalar<type, errorTrack>(n * cos(a), n * sin(a));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> ln(const ComplexScalar<type, errorTrack>& c) {
        const auto n = norm(c);
        const auto a = arg(c);
        return ComplexScalar<type, errorTrack>(ln(n), a);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> exp(const ComplexScalar<type, errorTrack>& c) {
        const auto& exp_real = exp(c.getReal());
        const auto& imagine = c.getImag();
        return ComplexScalar<type, errorTrack>(exp_real * cos(imagine), exp_real * sin(imagine));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> cos(const ComplexScalar<type, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<type, errorTrack>(cos(real) * cosh(imagine), - sin(real) * sinh(imagine));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sin(const ComplexScalar<type, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<type, errorTrack>(sin(real) * cosh(imagine), cos(real) * sinh(imagine));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> tan(const ComplexScalar<type, errorTrack>& c) {
        return sin(c) / cos(c);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sec(const ComplexScalar<type, errorTrack>& c) {
        return reciprocal(cos(c));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> csc(const ComplexScalar<type, errorTrack>& c) {
        return reciprocal(sin(c));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> cot(const ComplexScalar<type, errorTrack>& c) {
        return reciprocal(tan(c));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> cosh(const ComplexScalar<type, errorTrack>& c) {
        return (exp(c) + exp(-c)) >> 1;
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sinh(const ComplexScalar<type, errorTrack>& c) {
        return (exp(c) - exp(-c)) >> 1;
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> tanh(const ComplexScalar<type, errorTrack>& c) {
        return sinh(c) / cosh(c);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sech(const ComplexScalar<type, errorTrack>& c) {
        return reciprocal(cosh(c));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> csch(const ComplexScalar<type, errorTrack>& c) {
        return reciprocal(sinh(c));
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> coth(const ComplexScalar<type, errorTrack>& c) {
        return reciprocal(tanh(c));
    }
}

#endif