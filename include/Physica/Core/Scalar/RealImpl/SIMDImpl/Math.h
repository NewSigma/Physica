/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Scalar/ExprID.h"
#include "../SIMD.h"

namespace Physica {
    template<FloatPrec Prec, int Size>
    [[nodiscard]] __host__ __device__ auto unit(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>::select(x.isNegative(), SIMD<Real<Prec>, Size>(-1), SIMD<Real<Prec>, Size>(1));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] SIMD<Real<Prec>, Size> fma(const SIMD<Real<Prec>, Size> a, const SIMD<Real<Prec>, Size> b, const SIMD<Real<Prec>, Size> c) noexcept {
        return SIMD<Real<Prec>, Size>(mul_add(a.toMachine(), b.toMachine(), c.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] __host__ __device__ SIMD<Real<Prec>, Size> abs(SIMD<Real<Prec>, Size> x) noexcept {
        if constexpr (Prec == Float16)
            return SIMD<Real<Prec>, Size>(__habs2(x.toMachine()));
        else if constexpr (IsHost())
            return SIMD<Real<Prec>, Size>(abs(x.toMachine()));
        else
            noImpl();
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto relu(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>::select(x > SIMD<Real<Prec>, Size>::zeros(), x, SIMD<Real<Prec>, Size>::zeros());
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] SIMD<Real<Prec>, Size> square(SIMD<Real<Prec>, Size> x) noexcept {
        return x * x;
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto reciprocal(SIMD<Real<Prec>, Size> x) noexcept {
        // Vector partial reciprocal might divide by zero, do not assert it.
        using ResultType = SIMD<Real<Prec>, Size>;
        return ResultType(1) / x;
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] SIMD<Real<Prec>, Size> sqrt(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(sqrt(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] SIMD<Real<Prec>, Size> cbrt(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(cbrt(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto ln(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(log(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto ln1p(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(log1p(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto ln1pexp(SIMD<Real<Prec>, Size> x) noexcept {
        return relu(x) + ln1p(exp(-abs(x)));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto pow(SIMD<Real<Prec>, Size> x, SIMD<Real<Prec>, Size> y) noexcept {
        return SIMD<Real<Prec>, Size>(pow(x.toMachine(), y.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto sin(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(Physica::sin(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto cos(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(Physica::cos(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto sincos(SIMD<Real<Prec>, Size> x) noexcept {
        SIMD<Real<Prec>, Size> s, c;
        Physica::sincos(x.toMachine(), s.toMachine(), c.toMachine());
        return std::make_pair(s, c);
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto tan(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(Physica::tan(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto sec(SIMD<Real<Prec>, Size> x) noexcept {
        return reciprocal(cos(x));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto arctan2(SIMD<Real<Prec>, Size> y, SIMD<Real<Prec>, Size> x) noexcept {
        return Physica::atan2(y.toMachine(), x.toMachine());
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto tanh(SIMD<Real<Prec>, Size> x) noexcept {
        return SIMD<Real<Prec>, Size>(Physica::tanh(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] SIMD<Real<Prec>, Size> lncosh(SIMD<Real<Prec>, Size> x) noexcept {
        using Pack = SIMD<Real<Prec>, Size>;
        const auto x1 = abs(x);
        return x1 + ln1p(exp(-x1 * Pack(2))) - Pack(M_LN2);
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto round(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return Physica::round(x.toMachine());
    }
}

#include "Math/Exp.h"
