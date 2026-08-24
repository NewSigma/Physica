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
    [[nodiscard]] CoDiff<T> abs(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield abs(x_.value());
        auto& g = y.grad();
        x_.reverse(x_.isPositive() ? g : -g);
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> relu(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield relu(x_.value());
        x_.reverse(x_.isPositive() ? y.grad() : Tv(0));
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> square(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield square(x_.value());
        x_.reverse(Tv(2) * x_.value(), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> reciprocal(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield reciprocal(x_.value());
        x_.reverse(-square(y.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> sqrt(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield sqrt(x_.value());
        x_.reverse(Tv(0.5) / y.value(), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> cbrt(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield cbrt(x_.value());
        const auto x2_3 = y.value() / x_.value();
        x_.reverse(Tv(1.0 / 3) * x2_3, y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> ln(T&& x) noexcept requires(ReverseDiff<T>) {
        if constexpr (!T::isComplex())
            assert(x.isPositive() && "[Error]: Invalid param");

        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield ln(x_.value());
        x_.reverse(y.grad() / x_.value());
    }

    template<Scalar T>
    [[nodiscard]] auto ln1p(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield ln1p(x_.value());
        x_.reverse(y.grad() / (x_.value() + 1.0));
    }

    template<Scalar T>
    [[nodiscard]] auto log(const T& x, const T& a) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] CoDiff<T> exp(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield exp(x_.value());
        x_.reverse(y.value(), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> expm1(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield expm1(x_.value());
        x_.reverse(exp(x_.value()), y.grad());
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] CoDiff<T> pow(T&& x, U&& a) noexcept requires(ReverseDiff<T> && !Diffable<U>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        decltype(auto) a_ = decay_rvalue(std::forward<U>(a));
        auto& y = co_yield pow(x_.value(), a_);
        x_.reverse(y.value() * a_ / x_.value(), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] auto pow(const T& x, const T& n) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] CoDiff<T> cos(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        Tv c, s;
        sincos(x_.value(), s, c);
        auto& y = co_yield c;
        x_.reverse(-s, y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> sin(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto [s, c] = sincos(x_.value());
        auto& y = co_yield s;
        x_.reverse(c, y.grad());
    }

    template<Scalar T, Scalar U>
    [[nodiscard]] CoDiff<void> sincos(const T& x, U&& sin_result, U&& cos_result) noexcept requires(ReverseDiff<T> && ReverseDiff<U>) {
        using Tv = T::ValueType;
        Tv s, c;
        sincos(x.value(), s, c);
        decltype(auto) sin_ = decay_rvalue(std::forward<U>(sin_result));
        decltype(auto) cos_ = decay_rvalue(std::forward<U>(cos_result));
        sin_ = s;
        cos_ = c;
        co_await suspend_always{};
    }

    template<Scalar T>
    [[nodiscard]] auto tan(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto sec(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto csc(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto cot(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arccos(const T& x) noexcept requires(ReverseDiff<T>);

    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    [[nodiscard]] auto arcsin(const T& x) noexcept requires(ReverseDiff<T>);
    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    [[nodiscard]] auto arctan(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arcsec(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arccsc(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arccot(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] CoDiff<T> cosh(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield cosh(x_.value());
        x_.reverse(sinh(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> sinh(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield sinh(x_.value());
        x_.reverse(cosh(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> tanh(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield tanh(x_.value());
        x_.reverse(Tv(1) - square(y.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] auto sech(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto csch(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto coth(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arccosh(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arcsinh(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arctanh(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arcsech(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arccsch(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] auto arccoth(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] CoDiff<T> lncosh(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield lncosh(x_.value());
        x_.reverse(tanh(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> softplus(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield softplus(x_.value());
        x_.reverse(sigmoid(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] T floor(const T& x) noexcept requires(ReverseDiff<T>);

    template<Scalar T>
    [[nodiscard]] T ceil(const T& x) noexcept requires(ReverseDiff<T>);
}
