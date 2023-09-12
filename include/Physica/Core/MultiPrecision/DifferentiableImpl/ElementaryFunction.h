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
    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode> abs(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode> square(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode> reciprocal(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ Differentiable<ScalarType, Mode> sqrt(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> cbrt(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> ln(const Differentiable<ScalarType, Mode>& s);
/*
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> log(const Differentiable<ScalarType, Mode>& s, const Differentiable<ScalarType, Mode>& a);
*/
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> exp(const Differentiable<ScalarType, Mode>& s);
/*
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> pow(const Differentiable<ScalarType, Mode>& s, const Differentiable<ScalarType, Mode>& n);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> factorial(const Differentiable<ScalarType, Mode>& s);
*/
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> cos(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> sin(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    void sincos(Differentiable<ScalarType, Mode> s, Differentiable<ScalarType, Mode>& sin_result, Differentiable<ScalarType, Mode>& cos_result);
/*
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> tan(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> sec(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> csc(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> cot(const Differentiable<ScalarType, Mode>& s);
    //!Domain of definition: [0, Pi]
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arccos(const Differentiable<ScalarType, Mode>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arcsin(const Differentiable<ScalarType, Mode>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arctan(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arcsec(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arccsc(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arccot(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> cosh(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> sinh(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> tanh(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> sech(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> csch(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> coth(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arccosh(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arcsinh(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arctanh(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arcsech(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arccsch(const Differentiable<ScalarType, Mode>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode> arccoth(const Differentiable<ScalarType, Mode>& s);
*/
}

#include "ElementaryFunctionImpl.h"
