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
    __host__ __device__ inline Real<Option> abs(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> relu(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> square(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> reciprocal(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> sqrt(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> cbrt(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ln(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ln1p(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> log(const Real<Option>& s, const Real<Option>& a) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> exp(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> pow(const Real<Option>& s, const Real<Option>& n) noexcept;

    template<ScalarOption Option>
    Real<Option> factorial(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> cos(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> sin(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline void sincos(Real<Option> s, Real<Option>& sin_result, Real<Option>& cos_result) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> tan(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> sec(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> csc(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> cot(const Real<Option>& s) noexcept;
    //!Domain of definition: [0, Pi]
    template<ScalarOption Option>
    Real<Option> arccos(const Real<Option>& s) noexcept;
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption Option>
    Real<Option> arcsin(const Real<Option>& s) noexcept;
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<ScalarOption Option>
    Real<Option> arctan(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arcsec(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arccsc(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arccot(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> cosh(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> sinh(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> tanh(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> sech(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> csch(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> coth(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arccosh(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arcsinh(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arctanh(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arcsech(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arccsch(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> arccoth(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    Real<Option> lncosh(const Real<Option>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> floor(const Real<Option>& s) noexcept;
    
    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ceil(const Real<Option>& s) noexcept;
}

#include "MathImpl.h"
