/*
 * Copyright 2020-2024 Weibo He.
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

#include "../Complex.h"

namespace Physica {
    template<Scalar T>
    __host__ __device__ T abs(const Complex<T>& c);

    template<Scalar T>
    Complex<T> square(const Complex<T>& c);

    template<Scalar T>
    Complex<T> reciprocal(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sqrt(const Complex<T>& c);

    template<Scalar T>
    __host__ __device__ Complex<T> ln(const Complex<T>& c);

    template<Scalar T>
    __host__ __device__ Complex<T> ln1p(const Complex<T>& c);

    template<Scalar T>
    Complex<T> ln1pexp(const Complex<T>& c);

    template<Scalar T>
    Complex<T> exp(const Complex<T>& c);

    template<Scalar T>
    Complex<T> cos(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sin(const Complex<T>& c);

    template<Scalar T>
    Complex<T> tan(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sec(const Complex<T>& c);

    template<Scalar T>
    Complex<T> csc(const Complex<T>& c);

    template<Scalar T>
    Complex<T> cot(const Complex<T>& c);

    template<Scalar T>
    Complex<T> cosh(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sinh(const Complex<T>& c);

    template<Scalar T>
    Complex<T> tanh(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sech(const Complex<T>& c);

    template<Scalar T>
    Complex<T> csch(const Complex<T>& c);

    template<Scalar T>
    Complex<T> coth(const Complex<T>& c);

    template<Scalar T>
    __host__ __device__ Complex<T> lncosh(const Complex<T>& c);
}

#include "MathImpl.h"
