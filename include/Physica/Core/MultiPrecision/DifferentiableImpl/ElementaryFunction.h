/*
 * Copyright 2023 Weibo He.
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
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ inline Differentiable<ScalarType, Mode, Order> abs(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    inline Differentiable<ScalarType, Mode, Order> relu(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ inline Differentiable<ScalarType, Mode, Order> square(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ inline Differentiable<ScalarType, Mode, Order> reciprocal(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    __host__ __device__ Differentiable<ScalarType, Mode, Order> sqrt(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> cbrt(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> ln(const Differentiable<ScalarType, Mode, Order>& s);
/*
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> log(const Differentiable<ScalarType, Mode, Order>& s, const Differentiable<ScalarType, Mode, Order>& a);
*/
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> exp(const Differentiable<ScalarType, Mode, Order>& s);
/*
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> pow(const Differentiable<ScalarType, Mode, Order>& s, const Differentiable<ScalarType, Mode, Order>& n);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> factorial(const Differentiable<ScalarType, Mode, Order>& s);
*/
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> cos(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> sin(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    void sincos(Differentiable<ScalarType, Mode, Order> s, Differentiable<ScalarType, Mode, Order>& sin_result, Differentiable<ScalarType, Mode, Order>& cos_result);
/*
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> tan(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> sec(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> csc(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> cot(const Differentiable<ScalarType, Mode, Order>& s);
*/
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arccos(const Differentiable<ScalarType, Mode, Order>& s);
/*
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arcsin(const Differentiable<ScalarType, Mode, Order>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arctan(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arcsec(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arccsc(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arccot(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> cosh(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> sinh(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> tanh(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> sech(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> csch(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> coth(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arccosh(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arcsinh(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arctanh(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arcsech(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arccsch(const Differentiable<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, unsigned int Order>
    Differentiable<ScalarType, Mode, Order> arccoth(const Differentiable<ScalarType, Mode, Order>& s);
*/
}

#include "ElementaryFunctionImpl.h"
