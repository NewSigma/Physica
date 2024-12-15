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

#include "Physica/Core/MultiPrecision/Diff.h"

namespace Physica::Core {
    template<Scalar T>
    __host__ __device__ inline CoDiff<T> abs(T&& x) requires(Diffable<T>);

    template<Scalar T>
    inline CoDiff<T> relu(T&& x) requires(Diffable<T>);

    template<Scalar T>
    __host__ __device__ inline CoDiff<T> square(T&& x) requires(Diffable<T>);

    template<Scalar T>
    __host__ __device__ inline CoDiff<T> reciprocal(T&& x) requires(Diffable<T>);

    template<Scalar T>
    __host__ __device__ CoDiff<T> sqrt(T&& x) requires(Diffable<T>);

    template<Scalar T>
    CoDiff<T> cbrt(T&& x) requires(Diffable<T>);

    template<Scalar T>
    CoDiff<T> ln(T&& x) requires(Diffable<T>);

    template<Scalar T, int Order>
    auto ln1p(const Diff<T, DiffMode::Forward, Order>& x);
/*
    template<Scalar T, int Order>
    auto log(const Diff<T, DiffMode::Forward, Order>& x, const Diff<T, DiffMode::Forward, Order>& a);
*/
    template<Scalar T>
    CoDiff<T> exp(T&& x) requires(Diffable<T>);

    template<Scalar T, int Order>
    auto pow(const Diff<T, DiffMode::Forward, Order>& x, const T& a);
/*
    template<Scalar T, int Order>
    auto pow(const Diff<T, DiffMode::Forward, Order>& x, const Diff<T, DiffMode::Forward, Order>& n);
*/
    template<Scalar T>
    CoDiff<T> cos(T&& x) requires(Diffable<T>);

    template<Scalar T>
    CoDiff<T> sin(T&& x) requires(Diffable<T>);

    template<Scalar T, Scalar U>
    CoDiff<void> sincos(const T& x, U&& sin_result, U&& cos_result);

    template<Scalar T, int Order>
    auto tan(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto sec(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto csc(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto cot(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccos(const Diff<T, DiffMode::Forward, Order>& x);
/*
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T, int Order>
    auto arcsin(const Diff<T, DiffMode::Forward, Order>& x);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T, int Order>
    auto arctan(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arcsec(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccsc(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccot(const Diff<T, DiffMode::Forward, Order>& x);
*/
    template<Scalar T, int Order>
    auto cosh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto sinh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto tanh(const Diff<T, DiffMode::Forward, Order>& x);
/*
    template<Scalar T, int Order>
    auto sech(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto csch(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto coth(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccosh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arcsinh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arctanh(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arcsech(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccsch(const Diff<T, DiffMode::Forward, Order>& x);

    template<Scalar T, int Order>
    auto arccoth(const Diff<T, DiffMode::Forward, Order>& x);
*/
    template<Scalar T, int Order>
    auto lncosh(const Diff<T, DiffMode::Forward, Order>& x) noexcept;

    template<Scalar T, int Order>
    inline T floor(const Diff<T, DiffMode::Forward, Order>& x) {
        return floor(x.value());
    }
    
    template<Scalar T, int Order>
    inline T ceil(const Diff<T, DiffMode::Forward, Order>& x) {
        return ceil(x.value());
    }
}

#include "MathImpl.h"
