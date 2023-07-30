/*
 * Copyright 2019-2023 WeiBo He.
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
#pragma once
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<ScalarOption option>
    __host__ __device__ inline Scalar<option> abs(const Scalar<option>& s);

    template<ScalarOption option>
    __host__ __device__ Scalar<option> square(const Scalar<option>& s);

    template<ScalarOption option>
    __host__ __device__ inline Scalar<option> reciprocal(const Scalar<option>& s);

    template<ScalarOption option>
    __host__ __device__ Scalar<option> sqrt(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> cbrt(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> ln(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> log(const Scalar<option>& s, const Scalar<option>& a);

    template<ScalarOption option>
    Scalar<option> exp(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> pow(const Scalar<option>& s, const Scalar<option>& n);

    template<ScalarOption option>
    Scalar<option> factorial(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> cos(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> sin(const Scalar<option>& s);

    template<ScalarOption option>
    void sincos(Scalar<option> s, Scalar<option>& sin_result, Scalar<option>& cos_result);

    template<ScalarOption option>
    Scalar<option> tan(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> sec(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> csc(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> cot(const Scalar<option>& s);
    //!Domain of definition: [0, Pi]
    template<ScalarOption option>
    Scalar<option> arccos(const Scalar<option>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption option>
    Scalar<option> arcsin(const Scalar<option>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption option>
    Scalar<option> arctan(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arcsec(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arccsc(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arccot(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> cosh(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> sinh(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> tanh(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> sech(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> csch(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> coth(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arccosh(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arcsinh(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arctanh(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arcsech(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arccsch(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arccoth(const Scalar<option>& s);
}

#include "FunctionImpl/ElementaryImpl.h"
