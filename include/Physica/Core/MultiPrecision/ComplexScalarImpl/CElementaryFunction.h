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
#ifndef PHYSICA_CELEMENTARYFUNCTION_H
#define PHYSICA_CELEMENTARYFUNCTION_H
/*!
 * This file is part of implementations of \ComplexScalar.
 * Do not include this header file, include ComplexScalar.h instead.
 *
 * A number functions have not been implemented.
 */
namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> square(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    inline ComplexScalar<option, errorTrack> reciprocal(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sqrt(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> ln(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> exp(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> cos(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sin(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> tan(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sec(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> csc(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> cot(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> cosh(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sinh(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> tanh(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> sech(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> csch(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> coth(const ComplexScalar<option, errorTrack>& c);
}

#include "FunctionImpl/CElementaryImpl.h"

#endif