/*
 * Copyright 2019 WeiBo He.
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
#ifndef PHYSICA_ELEMENTARYFUNCTION_H
#define PHYSICA_ELEMENTARYFUNCTION_H
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> abs(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> square(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> reciprocal(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sqrt(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> cbrt(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> pow(const Scalar<option, errorTrack>& s, const Scalar<option, errorTrack>& n);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> factorial(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> ln(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    Scalar<option, errorTrack1 || errorTrack2> log(const Scalar<option, errorTrack1>& s, const Scalar<option, errorTrack2>& a);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> exp(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> cos(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sin(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> tan(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sec(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> csc(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> cot(const Scalar<option, errorTrack>& s);
    //!Domain of definition: [0, Pi]
    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccos(const Scalar<option, errorTrack>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arcsin(const Scalar<option, errorTrack>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arctan(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arcsec(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccsc(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccot(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> cosh(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sinh(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> tanh(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sech(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> csch(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> coth(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccosh(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arcsinh(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arctanh(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arcsech(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccsch(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccoth(const Scalar<option, errorTrack>& s);
}

#include "FunctionImpl/ElementaryImpl.h"

#endif