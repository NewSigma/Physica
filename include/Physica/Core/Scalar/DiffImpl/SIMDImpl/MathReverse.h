/*
 * Copyright 2024-2026 Weibo He.
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
    template<Packet T>
    [[nodiscard]] CoDiff<T> fma(T&& a, T&& b, T&& c) noexcept requires(ReverseDiff<T>) {
        decltype(auto) a_ = decay_rvalue(std::forward<T>(a));
        decltype(auto) b_ = decay_rvalue(std::forward<T>(b));
        decltype(auto) c_ = decay_rvalue(std::forward<T>(c));
        auto& y = co_yield fma(a_.value(), b_.value(), c_.value());
        auto& grad = y.grad();
        a_.reverse(b_.value(), grad);
        b_.reverse(a_.value(), grad);
        c_.reverse(grad);
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> abs(T&& x) noexcept requires(ReverseDiff<T>) {
        using Grad = std::remove_cvref_t<T>::GradType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield abs(x_.value());
        x_.reverse(Grad::select(x_.value().isPositive(), y.grad(), -y.grad()));
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> square(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_cvref_t<T>::ValueType::ScalarType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield square(x_.value());
        x_.reverse(Tv(2) * x_.value(), y.grad());
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> reciprocal(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield reciprocal(x_.value());
        x_.reverse(-square(y.value()), y.grad());
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> ln(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield ln(x_.value());
        x_.reverse(y.grad() / x_.value());
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> exp(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield exp(x_.value());
        x_.reverse(y.value(), y.grad());
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> expm1(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield expm1(x_.value());
        x_.reverse(exp(x_.value()), y.grad());
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> tanh(T&& x) noexcept requires(ReverseDiff<T>) {
        using Tv = std::remove_cvref_t<T>::ValueType::ScalarType;
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield tanh(x_.value());
        x_.reverse(Tv(1) - square(y.value()), y.grad());
    }

    template<Packet T>
    [[nodiscard]] CoDiff<T> lncosh(T&& x) noexcept requires(ReverseDiff<T>) {
        decltype(auto) x_ = decay_rvalue(std::forward<T>(x));
        auto& y = co_yield lncosh(x_.value());
        x_.reverse(tanh(x_.value()), y.grad());
    }
}
