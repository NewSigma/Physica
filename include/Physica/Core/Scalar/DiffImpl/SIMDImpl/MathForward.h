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
    [[nodiscard]] T fma(const T& a, const T& b, const T& c) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        auto value = fma(a.value(), b.value(), c.value());
        auto grad1 = fma(Grad(a), b.grad(), c.grad());
        auto grad2 = fma(Grad(b), a.grad(), grad1);
        return T(std::move(value), std::move(grad2));
    }

    template<Packet T>
    [[nodiscard]] T abs(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(abs(x.value()), Grad::select(x.value().isPositive(), x.grad(), -x.grad()));
    }

    template<Packet T>
    [[nodiscard]] T square(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(square(x.value()), Grad(x) * x.grad() * Grad(2));
    }

    template<Packet T>
    [[nodiscard]] T reciprocal(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const auto y = reciprocal(Grad(x));
        return T(y.value(), -x.grad() * square(y));
    }

    template<Packet T>
    [[nodiscard]] T ln(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(ln(x.value()), x.grad() / Grad(x));
    }

    template<Packet T>
    [[nodiscard]] T ln1p(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(ln1p(x.value()), x.grad() / (Grad(1) + Grad(x)));
    }

    template<Packet T>
    [[nodiscard]] T exp(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const auto y = exp(Grad(x));
        return T(y.value(), y * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T tanh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        const Grad y = tanh(Grad(x));
        return T(y.value(), (Grad(1) - square(y)) * x.grad());
    }

    template<Packet T>
    [[nodiscard]] T lncosh(const T& x) noexcept requires(ForwardDiff<T>) {
        using Grad = T::GradType;
        return T(lncosh(x.value()), tanh(Grad(x)) * x.grad());
    }
}
