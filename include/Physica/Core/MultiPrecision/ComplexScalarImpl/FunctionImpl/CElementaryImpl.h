/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PHYSICA_CELEMENTARYIMPL_H
#define PHYSICA_CELEMENTARYIMPL_H
/*!
 * This file is part of implementations of \ComplexScalar.
 * Do not include this header file, include ComplexScalar.h instead.
 */
namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> square(const ComplexScalar<option, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<option, errorTrack>(square(real) - square(imagine), (real * imagine) << 1);
    }

    template<ScalarOption option, bool errorTrack>
    inline ComplexScalar<option, errorTrack> reciprocal(const ComplexScalar<option, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        const auto divisor = reciprocal(square(real) + square(imagine));
        return ComplexScalar<option, errorTrack>(real * divisor, -imagine * divisor);
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sqrt(const ComplexScalar<option, errorTrack>& c) {
        const auto n = norm(c);
        const auto a = arg(c) >> 1;
        return ComplexScalar<option, errorTrack>(n * cos(a), n * sin(a));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> ln(const ComplexScalar<option, errorTrack>& c) {
        const auto n = norm(c);
        const auto a = arg(c);
        return ComplexScalar<option, errorTrack>(ln(n), a);
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> exp(const ComplexScalar<option, errorTrack>& c) {
        const auto& exp_real = exp(c.getReal());
        const auto& imagine = c.getImag();
        return ComplexScalar<option, errorTrack>(exp_real * cos(imagine), exp_real * sin(imagine));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> cos(const ComplexScalar<option, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<option, errorTrack>(cos(real) * cosh(imagine), - sin(real) * sinh(imagine));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sin(const ComplexScalar<option, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImag();
        return ComplexScalar<option, errorTrack>(sin(real) * cosh(imagine), cos(real) * sinh(imagine));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> tan(const ComplexScalar<option, errorTrack>& c) {
        return sin(c) / cos(c);
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sec(const ComplexScalar<option, errorTrack>& c) {
        return reciprocal(cos(c));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> csc(const ComplexScalar<option, errorTrack>& c) {
        return reciprocal(sin(c));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> cot(const ComplexScalar<option, errorTrack>& c) {
        return reciprocal(tan(c));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> cosh(const ComplexScalar<option, errorTrack>& c) {
        return (exp(c) + exp(-c)) >> 1;
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sinh(const ComplexScalar<option, errorTrack>& c) {
        return (exp(c) - exp(-c)) >> 1;
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> tanh(const ComplexScalar<option, errorTrack>& c) {
        return sinh(c) / cosh(c);
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sech(const ComplexScalar<option, errorTrack>& c) {
        return reciprocal(cosh(c));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> csch(const ComplexScalar<option, errorTrack>& c) {
        return reciprocal(sinh(c));
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> coth(const ComplexScalar<option, errorTrack>& c) {
        return reciprocal(tanh(c));
    }
}

#endif