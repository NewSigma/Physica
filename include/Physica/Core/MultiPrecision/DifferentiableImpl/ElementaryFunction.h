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
    __host__ __device__ inline Differentiable<ScalarType, Mode, 1> abs(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    inline Differentiable<ScalarType, Mode, 1> relu(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode, 1> square(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ inline Differentiable<ScalarType, Mode, 1> reciprocal(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    __host__ __device__ Differentiable<ScalarType, Mode, 1> sqrt(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> cbrt(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> ln(const Differentiable<ScalarType, Mode, 1>& s);
/*
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> log(const Differentiable<ScalarType, Mode, 1>& s, const Differentiable<ScalarType, Mode, 1>& a);
*/
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> exp(const Differentiable<ScalarType, Mode, 1>& s);
/*
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> pow(const Differentiable<ScalarType, Mode, 1>& s, const Differentiable<ScalarType, Mode, 1>& n);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> factorial(const Differentiable<ScalarType, Mode, 1>& s);
*/
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> cos(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> sin(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    void sincos(Differentiable<ScalarType, Mode, 1> s, Differentiable<ScalarType, Mode, 1>& sin_result, Differentiable<ScalarType, Mode, 1>& cos_result);
/*
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> tan(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> sec(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> csc(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> cot(const Differentiable<ScalarType, Mode, 1>& s);
*/
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arccos(const Differentiable<ScalarType, Mode, 1>& s);
/*
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arcsin(const Differentiable<ScalarType, Mode, 1>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arctan(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arcsec(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arccsc(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arccot(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> cosh(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> sinh(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> tanh(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> sech(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> csch(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> coth(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arccosh(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arcsinh(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arctanh(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arcsech(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arccsch(const Differentiable<ScalarType, Mode, 1>& s);

    template<class ScalarType, DiffMode Mode>
    Differentiable<ScalarType, Mode, 1> arccoth(const Differentiable<ScalarType, Mode, 1>& s);
*/
}

#include "ElementaryFunctionImpl.h"
