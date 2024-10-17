/*
 * Copyright 2019-2024 Weibo He.
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
    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> abs(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> relu(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> square(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> reciprocal(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> sqrt(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> cbrt(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ln(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ln1p(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> log(const Scalar<Option>& s, const Scalar<Option>& a) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> exp(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> pow(const Scalar<Option>& s, const Scalar<Option>& n) noexcept;

    template<ScalarOption Option>
    Scalar<Option> factorial(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> cos(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> sin(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline void sincos(Scalar<Option> s, Scalar<Option>& sin_result, Scalar<Option>& cos_result) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> tan(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> sec(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> csc(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> cot(const Scalar<Option>& s) noexcept;
    //!Domain of definition: [0, Pi]
    template<ScalarOption Option>
    Scalar<Option> arccos(const Scalar<Option>& s) noexcept;
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption Option>
    Scalar<Option> arcsin(const Scalar<Option>& s) noexcept;
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption Option>
    Scalar<Option> arctan(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arcsec(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arccsc(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arccot(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> cosh(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> sinh(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> tanh(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> sech(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> csch(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> coth(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arccosh(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arcsinh(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arctanh(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arcsech(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arccsch(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arccoth(const Scalar<Option>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> lncosh(const Scalar<Option>& s) noexcept;
}

#include "FunctionImpl/ElementaryImpl.h"
