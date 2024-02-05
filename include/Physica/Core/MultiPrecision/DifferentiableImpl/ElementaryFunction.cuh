/*
 * Copyright 2024 WeiBo He.
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
    template<class ScalarType>
    inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>> abs(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>> relu(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>> square(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>> reciprocal(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> sqrt(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> cbrt(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);
*/
    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> ln(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);
/*
    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> log(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s, const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& a);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> exp(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> pow(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s, const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& n);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> factorial(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> cos(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> sin(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    void sincos(device_obj<Differentiable<ScalarType, DiffMode::Reverse>> s, device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& sin_result, device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& cos_result);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> tan(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> sec(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> csc(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> cot(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arccos(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arcsin(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arctan(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arcsec(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arccsc(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arccot(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> cosh(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> sinh(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> tanh(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> sech(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> csch(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> coth(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arccosh(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arcsinh(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arctanh(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arcsech(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arccsch(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>> arccoth(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s);
*/
}

#include "ElementaryFunctionImpl.cuh"
