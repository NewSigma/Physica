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
    template<Scalar T, int Size>
    [[nodiscard]] __host__ __device__ auto unit(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>::select(x.isNegative(), SIMD<T, Size>(-1), SIMD<T, Size>(1));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> fma(const SIMD<T, Size> a, const SIMD<T, Size> b, const SIMD<T, Size> c) noexcept {
        return SIMD<T, Size>(mul_add(a.toMachine(), b.toMachine(), c.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] __host__ __device__ SIMD<T, Size> abs(SIMD<T, Size> x) noexcept {
        if constexpr (T::Prec == Float16)
            return SIMD<T, Size>(__habs2(x.toMachine()));
        else if constexpr (IsHost())
            return SIMD<T, Size>(abs(x.toMachine()));
        else
            noImpl();
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto relu(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>::select(x > SIMD<T, Size>::zeros(), x, SIMD<T, Size>::zeros());
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> square(SIMD<T, Size> x) noexcept {
        return x * x;
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto reciprocal(SIMD<T, Size> x) noexcept {
        // Vector partial reciprocal might divide by zero, do not assert it.
        using ResultType = SIMD<T, Size>;
        return ResultType(1) / x;
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> sqrt(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(sqrt(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> cbrt(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(cbrt(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto ln(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(log(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto ln1p(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(log1p(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto ln1pexp(SIMD<T, Size> x) noexcept {
        return relu(x) + ln1p(exp(-abs(x)));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto pow(SIMD<T, Size> x, SIMD<T, Size> y) noexcept {
        return SIMD<T, Size>(pow(x.toMachine(), y.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto sin(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(Physica::sin(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto cos(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(Physica::cos(x.toMachine()));
    }

    template<Scalar T, int Size>
    void sincos(SIMD<T, Size> x, SIMD<T, Size>& s, SIMD<T, Size>& c) noexcept {
        Physica::sincos(x.toMachine(), s.toMachine(), c.toMachine());
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto tan(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(Physica::tan(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto sec(SIMD<T, Size> x) noexcept {
        return reciprocal(cos(x));
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto arctan2(SIMD<T, Size> y, SIMD<T, Size> x) noexcept {
        return Physica::atan2(y.toMachine(), x.toMachine());
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto tanh(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(Physica::tanh(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> lncosh(SIMD<T, Size> x) noexcept {
        using Pack = SIMD<T, Size>;
        const auto x1 = abs(x);
        return x1 + ln1p(exp(-x1 * Pack(2))) - Pack(M_LN2);
    }

    template<Scalar T, int Size>
    [[nodiscard]] auto round(SIMD<T, Size> x) noexcept -> SIMD<T, Size> {
        return Physica::round(x.toMachine());
    }
}

#include "Math/Exp.h"
