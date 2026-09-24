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

    template<Scalar T, Scalar U>
    [[nodiscard]] CoDiff<T> log(T&& x, U&& a) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        decltype(auto) a_ = decay_rvalue(std::forward<U>(a));
        auto& y = co_yield log(x_.value(), a_.value());
        const auto lna = ln(a_.value());
        const auto& g = y.grad();
        x_.reverse(reciprocal(x_.value() * lna), g);
        if constexpr (ReverseDiff<U>)
            a_.reverse(-ln(x_.value()) / (a_.value() * square(lna)), g);
    }

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

    template<Scalar T, Scalar U>
    [[nodiscard]] CoDiff<T> pow(T&& x, U&& n) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        decltype(auto) n_ = decay_rvalue(std::forward<U>(n));
        auto& y = co_yield pow(x_.value(), n_.value());
        const auto& g = y.grad();
        x_.reverse(y.value() * n_.value() / x_.value(), g);
        if constexpr (ReverseDiff<U>)
            n_.reverse(y.value() * ln(x_.value()), g);
    }

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
    [[nodiscard]] CoDiff<T> tan(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield tan(x_.value());
        x_.reverse(square(sec(x_.value())), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> sec(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield sec(x_.value());
        x_.reverse(sec(x_.value()) * tan(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> csc(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield csc(x_.value());
        x_.reverse(-csc(x_.value()) * cot(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> cot(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield cot(x_.value());
        x_.reverse(-square(csc(x_.value())), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arccos(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arccos(x_.value());
        x_.reverse(-reciprocal(sqrt(Tv(1) - square(x_.value()))), y.grad());
    }

    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    [[nodiscard]] CoDiff<T> arcsin(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arcsin(x_.value());
        x_.reverse(reciprocal(sqrt(Tv(1) - square(x_.value()))), y.grad());
    }

    //!Domain of definition: [-Pi / 2, Pi / 2]
    template<Scalar T>
    [[nodiscard]] CoDiff<T> arctan(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arctan(x_.value());
        x_.reverse(reciprocal(Tv(1) + square(x_.value())), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arcsec(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arcsec(x_.value());
        const auto u = reciprocal(x_.value());
        x_.reverse(reciprocal(square(x_.value()) * sqrt(Tv(1) - square(u))), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arccsc(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arccsc(x_.value());
        const auto u = reciprocal(x_.value());
        x_.reverse(-reciprocal(square(x_.value()) * sqrt(Tv(1) - square(u))), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arccot(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arccot(x_.value());
        x_.reverse(-reciprocal(Tv(1) + square(x_.value())), y.grad());
    }

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
    [[nodiscard]] CoDiff<T> sech(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield sech(x_.value());
        x_.reverse(-sech(x_.value()) * tanh(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> csch(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield csch(x_.value());
        x_.reverse(-csch(x_.value()) * coth(x_.value()), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> coth(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield coth(x_.value());
        x_.reverse(-square(csch(x_.value())), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arccosh(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arccosh(x_.value());
        x_.reverse(reciprocal(sqrt(square(x_.value()) - Tv(1))), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arcsinh(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arcsinh(x_.value());
        x_.reverse(reciprocal(sqrt(square(x_.value()) + Tv(1))), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arctanh(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arctanh(x_.value());
        x_.reverse(reciprocal(Tv(1) - square(x_.value())), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arcsech(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arcsech(x_.value());
        const auto u = reciprocal(x_.value());
        x_.reverse(-reciprocal(square(x_.value()) * sqrt(square(u) - Tv(1))), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arccsch(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arccsch(x_.value());
        const auto u = reciprocal(x_.value());
        x_.reverse(-reciprocal(square(x_.value()) * sqrt(square(u) + Tv(1))), y.grad());
    }

    template<Scalar T>
    [[nodiscard]] CoDiff<T> arccoth(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_reference_t<T>::ValueType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield arccoth(x_.value());
        x_.reverse(reciprocal(Tv(1) - square(x_.value())), y.grad());
    }

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
    [[nodiscard]] T floor(const T& x) noexcept requires(ReverseDiff<T>) {
        return T(floor(x.value()));
    }

    template<Scalar T>
    [[nodiscard]] T ceil(const T& x) noexcept requires(ReverseDiff<T>) {
        return T(ceil(x.value()));
    }
}
