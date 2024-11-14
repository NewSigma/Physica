/*
 * Copyright 2023-2024 Weibo He.
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
    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> abs(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    inline Diff<ScalarType, Mode, Order> relu(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> square(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<ScalarType, Mode, Order> reciprocal(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    __host__ __device__ Diff<ScalarType, Mode, Order> sqrt(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cbrt(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> ln(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> ln1p(const Diff<ScalarType, Mode, Order>& x);
/*
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> log(const Diff<ScalarType, Mode, Order>& s, const Diff<ScalarType, Mode, Order>& a);
*/
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> exp(const Diff<ScalarType, Mode, Order>& s);
/*
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> pow(const Diff<ScalarType, Mode, Order>& s, const Diff<ScalarType, Mode, Order>& n);
*/
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cos(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> sin(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    void sincos(Diff<ScalarType, Mode, Order> s, Diff<ScalarType, Mode, Order>& sin_result, Diff<ScalarType, Mode, Order>& cos_result);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> tan(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> sec(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> csc(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cot(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arccos(const Diff<ScalarType, Mode, Order>& s);
/*
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arcsin(const Diff<ScalarType, Mode, Order>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arctan(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arcsec(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arccsc(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arccot(const Diff<ScalarType, Mode, Order>& s);
*/
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> cosh(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> sinh(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> tanh(const Diff<ScalarType, Mode, Order>& s);
/*
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> sech(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> csch(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> coth(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arccosh(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arcsinh(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arctanh(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arcsech(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arccsch(const Diff<ScalarType, Mode, Order>& s);

    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> arccoth(const Diff<ScalarType, Mode, Order>& s);
*/
    template<class ScalarType, DiffMode Mode, int Order>
    Diff<ScalarType, Mode, Order> lncosh(const Diff<ScalarType, Mode, Order>& s) noexcept;

    template<class ScalarType, DiffMode Mode, int Order>
    inline ScalarType floor(const Diff<ScalarType, Mode, Order>& s) {
        return floor(s.getValue());
    }
    
    template<class ScalarType, DiffMode Mode, int Order>
    inline ScalarType ceil(const Diff<ScalarType, Mode, Order>& s) {
        return ceil(s.getValue());
    }
}

#include "MathImpl.h"
