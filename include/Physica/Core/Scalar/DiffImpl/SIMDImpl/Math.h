/*
 * Copyright 2024-2025 Weibo He.
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
    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] CoDiff<SIMD<Diff<T, Mode, Order>, Size>> fma(
            const SIMD<Diff<T, Mode, Order>, Size>& a,
            const SIMD<Diff<T, Mode, Order>, Size>& b,
            const SIMD<Diff<T, Mode, Order>, Size>& c) noexcept {
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradType = ResultType::GradType;
            auto value = fma(a.value(), b.value(), c.value());
            auto grad1 = fma(GradType(a), b.grad(), c.grad());
            auto grad2 = fma(GradType(b), a.grad(), grad1);
            co_return ResultType(std::move(value), std::move(grad2));
        }
        else {
            auto result = ResultType(fma(a.value(), b.value(), c.value()));
            co_yield result;

            auto& grad = result.grad();
            a.reverse(grad * b.value());
            b.reverse(grad * a.value());
            c.reverse(grad);
        }
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto abs(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(abs(x.value()), GradPacket::select(x.value().isPositive(), x.grad(), -x.grad()));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto square(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(square(x.value()), GradPacket(x) * x.grad() * T(2));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto reciprocal(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        const auto y = reciprocal(GradPacket(x));
        return ResultType(y.value(), -x.grad() * square(y));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto ln(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(ln(x.value()), x.grad() / GradPacket(x));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto ln1p(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(ln1p(x.value()), x.grad() / (GradPacket(1) + GradPacket(x)));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto exp(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        const auto y = exp(GradPacket(x));
        return ResultType(y.value(), y * x.grad());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto tanh(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        const GradPacket y = tanh(GradPacket(x));
        return ResultType(y.value(), (GradPacket(1) - square(y)) * x.grad());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    [[nodiscard]] auto lncosh(const SIMD<Diff<T, Mode, Order>, Size>& x) noexcept {
        static_assert(Mode != DiffMode::Reverse, "[Error]: Not implemented");
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        using GradPacket = ResultType::GradType;
        return ResultType(lncosh(x.value()), tanh(GradPacket(x)) * x.grad());
    }
}
