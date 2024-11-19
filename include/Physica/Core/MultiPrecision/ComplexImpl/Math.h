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

namespace Physica::Core {
    template<Scalar T>
    T abs(const Complex<T>& c);

    template<Scalar T>
    Complex<T> square(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> reciprocal(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sqrt(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> ln(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> exp(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> cos(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> sin(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> tan(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sec(const Complex<T>& c);

    template<Scalar T>
    Complex<T> csc(const Complex<T>& c);

    template<Scalar T>
    Complex<T> cot(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> cosh(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> sinh(const Complex<T>& c);

    template<Scalar T>
    inline Complex<T> tanh(const Complex<T>& c);

    template<Scalar T>
    Complex<T> sech(const Complex<T>& c);

    template<Scalar T>
    Complex<T> csch(const Complex<T>& c);

    template<Scalar T>
    Complex<T> coth(const Complex<T>& c);

    template<Scalar T>
    Complex<T> lncosh(const Complex<T>& c);
}

#include "MathImpl.h"
