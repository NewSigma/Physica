/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Scalar/Diff.h"

namespace Physica {
    template<Scalar T>
    inline CoDiff<T> abs(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield abs(x_.value());
        auto& g = y.grad();
        x_.reverse(x_.isPositive() ? g : -g);
    }

    template<Scalar T>
    inline CoDiff<T> relu(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield relu(x_.value());
        x_.reverse(x_.isPositive() ? y.grad() : Tv(0));
    }

    template<Scalar T>
    inline CoDiff<T> square(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield square(x_.value());
        x_.reverse(Tv(2) * x_.value() * y.grad());
    }

    template<Scalar T>
    inline CoDiff<T> reciprocal(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield reciprocal(x_.value());
        x_.reverse(-square(y.value()) * y.grad());
    }

    template<Scalar T>
    CoDiff<T> sqrt(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield sqrt(x_.value());
        x_.reverse(y.grad() / y.value() * Tv(0.5));
    }

    template<Scalar T>
    CoDiff<T> cbrt(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield cbrt(x_.value());
        const auto x2_3 = y.value() / x_.value();
        x_.reverse((Tv(1.0 / 3) * y.grad()) * x2_3);
    }

    template<Scalar T>
    CoDiff<T> ln(T&& x) noexcept requires(ReverseDiff<T>) {
        if constexpr (!T::isComplex)
            assert(x.isPositive() && "[Error]: Invalid param");

        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield ln(x_.value());
        x_.reverse(y.grad() / x_.value());
    }

    // ln1p
    // ln1pexp

    //template<Scalar T, int Order>
    //auto log(const Diff<T, DiffMode::Forward, Order>& x, const Diff<T, DiffMode::Forward, Order>& a);
    
    template<Scalar T>
    CoDiff<T> exp(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield exp(x_.value());
        x_.reverse(y.value() * y.grad());
    }

    template<Scalar T, Scalar U>
    CoDiff<T> pow(T&& x, U&& a) noexcept requires(ReverseDiff<T> && !Diffable<U>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        LazyDestroy<U&&> a_ = std::forward<U>(a);
        auto y = co_yield pow(x_.value(), a_);
        x_.reverse(y.value() * y.grad() * a_ / x_.value());
    }

    //template<Scalar T, int Order>
    //auto pow(const Diff<T, DiffMode::Forward, Order>& x, const Diff<T, DiffMode::Forward, Order>& n);

    template<Scalar T>
    CoDiff<T> cos(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        Tv c, s;
        sincos(x_.value(), s, c);
        auto y = co_yield c;
        x_.reverse(-s * y.grad());
    }

    template<Scalar T>
    CoDiff<T> sin(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        Tv c, s;
        sincos(x_.value(), s, c);
        auto y = co_yield s;
        x_.reverse(c * y.grad());
    }

    template<Scalar T, Scalar U>
    CoDiff<void> sincos(const T& x, U&& sin_result, U&& cos_result) noexcept requires(ReverseDiff<T> && ReverseDiff<U>) {
        using Tv = T::ValueType;
        Tv s, c;
        sincos(x.value(), s, c);
        LazyDestroy<U&&> sin_ = std::forward<U>(sin_result);
        LazyDestroy<U&&> cos_ = std::forward<U>(cos_result);
        sin_ = s;
        cos_ = c;
        co_await std::suspend_always{};
    }
    /*
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

    template<Scalar T>
    CoDiff<T> cosh(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield cosh(x_.value());
        x_.reverse(sinh(x_.value()) * y.grad());
    }

    template<Scalar T>
    CoDiff<T> sinh(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield sinh(x_.value());
        x_.reverse(cosh(x_.value()) * y.grad());
    }

    template<Scalar T>
    CoDiff<T> tanh(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield tanh(x_.value());
        x_.reverse((Tv(1) - square(y.value())) * y.grad());
    }

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

    template<Scalar T>
    CoDiff<T> lncosh(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto y = co_yield lncosh(x_.value());
        x_.reverse(tanh(x_.value()) * y.grad());
    }

    /*
    template<Scalar T, int Order>
    inline T floor(const Diff<T, DiffMode::Forward, Order>& x);
    
    template<Scalar T, int Order>
    inline T ceil(const Diff<T, DiffMode::Forward, Order>& x);
    */
}
