/*
 * Copyright 2025-2026 Weibo He.
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
    CoDiff<T> abs(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield abs(x_.value());
        auto& g = y.grad();
        x_.reverse(x_.isPositive() ? g : -g);
    }

    template<Scalar T>
    CoDiff<T> relu(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield relu(x_.value());
        x_.reverse(x_.isPositive() ? y.grad() : Tv(0));
    }

    template<Scalar T>
    CoDiff<T> square(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield square(x_.value());
        x_.reverse(Tv(2) * x_.value() * y.grad());
    }

    template<Scalar T>
    CoDiff<T> reciprocal(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield reciprocal(x_.value());
        x_.reverse(-square(y.value()) * y.grad());
    }

    template<Scalar T>
    CoDiff<T> sqrt(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield sqrt(x_.value());
        x_.reverse(y.grad() / y.value() * Tv(0.5));
    }

    template<Scalar T>
    CoDiff<T> cbrt(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield cbrt(x_.value());
        const auto x2_3 = y.value() / x_.value();
        x_.reverse((Tv(1.0 / 3) * y.grad()) * x2_3);
    }

    template<Scalar T>
    CoDiff<T> ln(T&& x) noexcept requires(ReverseDiff<T>) {
        if constexpr (!T::isComplex)
            assert(x.isPositive() && "[Error]: Invalid param");

        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield ln(x_.value());
        x_.reverse(y.grad() / x_.value());
    }

    template<Scalar T>
    auto ln1p(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto ln1pexp(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto log(const T& x, const T& a) noexcept requires(ReverseDiff<T>);
    
    template<Scalar T>
    CoDiff<T> exp(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield exp(x_.value());
        x_.reverse(y.value() * y.grad());
    }

    template<Scalar T, Scalar U>
    CoDiff<T> pow(T&& x, U&& a) noexcept requires(ReverseDiff<T> && !Diffable<U>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        LazyDestroy<U&&> a_ = std::forward<U>(a);
        auto& y = co_yield pow(x_.value(), a_);
        x_.reverse(y.value() * y.grad() * a_ / x_.value());
    }

    template<Scalar T>
    auto pow(const T& x, const T& n) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    CoDiff<T> cos(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        Tv c, s;
        sincos(x_.value(), s, c);
        auto& y = co_yield c;
        x_.reverse(-s * y.grad());
    }

    template<Scalar T>
    CoDiff<T> sin(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        Tv c, s;
        sincos(x_.value(), s, c);
        auto& y = co_yield s;
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
        co_await suspend_always{};
    }

    template<Scalar T>
    auto tan(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto sec(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto csc(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto cot(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arccos(const T& x) noexcept requires(ReverseDiff<T>);

    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    auto arcsin(const T& x) noexcept requires(ReverseDiff<T>);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    auto arctan(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arcsec(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arccsc(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arccot(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    CoDiff<T> cosh(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield cosh(x_.value());
        x_.reverse(sinh(x_.value()) * y.grad());
    }

    template<Scalar T>
    CoDiff<T> sinh(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield sinh(x_.value());
        x_.reverse(cosh(x_.value()) * y.grad());
    }

    template<Scalar T>
    CoDiff<T> tanh(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield tanh(x_.value());
        x_.reverse((Tv(1) - square(y.value())) * y.grad());
    }

    template<Scalar T>
    auto sech(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto csch(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto coth(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arccosh(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arcsinh(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arctanh(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arcsech(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arccsch(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    auto arccoth(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    CoDiff<T> lncosh(T&& x) noexcept requires(ReverseDiff<T>) {
        LazyDestroy<T&&> x_ = std::forward<T>(x);
        auto& y = co_yield lncosh(x_.value());
        x_.reverse(tanh(x_.value()) * y.grad());
    }

    template<Scalar T>
    T floor(const T& x) noexcept requires(ReverseDiff<T>);
    
    template<Scalar T>
    T ceil(const T& x) noexcept requires(ReverseDiff<T>);
}
