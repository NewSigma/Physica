/*
 * Copyright 2024 Weibo He.
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
/*
    template<Scalar T>
    inline device_obj<Diff<T, DiffMode::Reverse, 1>> abs(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    inline device_obj<Diff<T, DiffMode::Reverse, 1>> relu(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    inline device_obj<Diff<T, DiffMode::Reverse, 1>> square(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    inline device_obj<Diff<T, DiffMode::Reverse, 1>> reciprocal(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> sqrt(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> cbrt(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);
*/
    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> ln(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);
/*
    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> log(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s, const device_obj<Diff<T, DiffMode::Reverse, 1>>& a);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> exp(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> pow(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s, const device_obj<Diff<T, DiffMode::Reverse, 1>>& n);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> factorial(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> cos(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> sin(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    void sincos(device_obj<Diff<T, DiffMode::Reverse, 1>> s, device_obj<Diff<T, DiffMode::Reverse, 1>>& sin_result, device_obj<Diff<T, DiffMode::Reverse, 1>>& cos_result);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> tan(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> sec(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> csc(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> cot(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arccos(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arcsin(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arctan(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arcsec(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arccsc(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arccot(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> cosh(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> sinh(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> tanh(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> sech(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> csch(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> coth(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arccosh(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arcsinh(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arctanh(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arcsech(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arccsch(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);

    template<Scalar T>
    device_obj<Diff<T, DiffMode::Reverse, 1>> arccoth(const device_obj<Diff<T, DiffMode::Reverse, 1>>& s);
*/
}

#include "MathImpl.cuh"
