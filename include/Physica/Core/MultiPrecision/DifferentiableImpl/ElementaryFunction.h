/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType>
    __host__ __device__ inline Differentiable<ScalarType> abs(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    __host__ __device__ inline Differentiable<ScalarType> square(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    __host__ __device__ inline Differentiable<ScalarType> reciprocal(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    __host__ __device__ Differentiable<ScalarType> sqrt(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> cbrt(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> ln(const Differentiable<ScalarType>& s);
/*
    template<class ScalarType>
    Differentiable<ScalarType> log(const Differentiable<ScalarType>& s, const Differentiable<ScalarType>& a);
*/
    template<class ScalarType>
    Differentiable<ScalarType> exp(const Differentiable<ScalarType>& s);
/*
    template<class ScalarType>
    Differentiable<ScalarType> pow(const Differentiable<ScalarType>& s, const Differentiable<ScalarType>& n);

    template<class ScalarType>
    Differentiable<ScalarType> factorial(const Differentiable<ScalarType>& s);
*/
    template<class ScalarType>
    Differentiable<ScalarType> cos(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> sin(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    void sincos(Differentiable<ScalarType> s, Differentiable<ScalarType>& sin_result, Differentiable<ScalarType>& cos_result);
/*
    template<class ScalarType>
    Differentiable<ScalarType> tan(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> sec(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> csc(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> cot(const Differentiable<ScalarType>& s);
    //!Domain of definition: [0, Pi]
    template<class ScalarType>
    Differentiable<ScalarType> arccos(const Differentiable<ScalarType>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType>
    Differentiable<ScalarType> arcsin(const Differentiable<ScalarType>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType>
    Differentiable<ScalarType> arctan(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arcsec(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arccsc(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arccot(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> cosh(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> sinh(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> tanh(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> sech(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> csch(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> coth(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arccosh(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arcsinh(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arctanh(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arcsech(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arccsch(const Differentiable<ScalarType>& s);

    template<class ScalarType>
    Differentiable<ScalarType> arccoth(const Differentiable<ScalarType>& s);
*/
}

#include "ElementaryFunctionImpl.h"
