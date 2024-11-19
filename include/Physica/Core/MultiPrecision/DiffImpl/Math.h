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
    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<T, Mode, Order> abs(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    inline Diff<T, Mode, Order> relu(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<T, Mode, Order> square(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline Diff<T, Mode, Order> reciprocal(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ Diff<T, Mode, Order> sqrt(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> cbrt(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> ln(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> ln1p(const Diff<T, Mode, Order>& x);
/*
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> log(const Diff<T, Mode, Order>& s, const Diff<T, Mode, Order>& a);
*/
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> exp(const Diff<T, Mode, Order>& s);
/*
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> pow(const Diff<T, Mode, Order>& s, const Diff<T, Mode, Order>& n);
*/
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> cos(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> sin(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    void sincos(Diff<T, Mode, Order> s, Diff<T, Mode, Order>& sin_result, Diff<T, Mode, Order>& cos_result);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> tan(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> sec(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> csc(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> cot(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arccos(const Diff<T, Mode, Order>& s);
/*
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arcsin(const Diff<T, Mode, Order>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arctan(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arcsec(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arccsc(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arccot(const Diff<T, Mode, Order>& s);
*/
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> cosh(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> sinh(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> tanh(const Diff<T, Mode, Order>& s);
/*
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> sech(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> csch(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> coth(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arccosh(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arcsinh(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arctanh(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arcsech(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arccsch(const Diff<T, Mode, Order>& s);

    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> arccoth(const Diff<T, Mode, Order>& s);
*/
    template<Scalar T, DiffMode Mode, int Order>
    Diff<T, Mode, Order> lncosh(const Diff<T, Mode, Order>& s) noexcept;

    template<Scalar T, DiffMode Mode, int Order>
    inline T floor(const Diff<T, Mode, Order>& s) {
        return floor(s.getValue());
    }
    
    template<Scalar T, DiffMode Mode, int Order>
    inline T ceil(const Diff<T, Mode, Order>& s) {
        return ceil(s.getValue());
    }
}

#include "MathImpl.h"
