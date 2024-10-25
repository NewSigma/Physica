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
    template<class ScalarType>
    ScalarType abs(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> square(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> reciprocal(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> sqrt(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> ln(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> exp(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> cos(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> sin(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> tan(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> sec(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> csc(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> cot(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> cosh(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> sinh(const Complex<ScalarType>& c);

    template<class ScalarType>
    inline Complex<ScalarType> tanh(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> sech(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> csch(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> coth(const Complex<ScalarType>& c);

    template<class ScalarType>
    Complex<ScalarType> lncosh(const Complex<ScalarType>& c);
}

#include "MathImpl.h"
