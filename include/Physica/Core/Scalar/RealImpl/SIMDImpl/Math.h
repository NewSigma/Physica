/*
 * Copyright 2023-2025 Weibo He.
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
    [[nodiscard]] __host__ __device__ SIMD<T, Size> abs(SIMD<T, Size> x) noexcept {
        if constexpr (T::Prec == Float16)
            return SIMD<T, Size>(__habs2(x.toMachine()));
        else if constexpr (IsHost())
            return SIMD<T, Size>(abs(x.toMachine()));
        else
            noImpl(__func__);
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
    auto cos(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(Physica::cos(x.toMachine()));
    }

    template<Scalar T, int Size>
    void sincos(SIMD<T, Size> x, SIMD<T, Size>& s, SIMD<T, Size>& c) noexcept {
        Physica::sincos(x.toMachine(), s.toMachine(), c.toMachine());
    }

    template<Scalar T, int Size>
    auto tan(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(Physica::tan(x.toMachine()));
    }

    template<Scalar T, int Size>
    auto sec(SIMD<T, Size> x) noexcept {
        return reciprocal(cos(x));
    }

    template<Scalar T, int Size>
    auto arctan2(SIMD<T, Size> y, SIMD<T, Size> x) noexcept {
        return Physica::atan2(y.toMachine(), x.toMachine());
    }

    template<Scalar T, int Size>
    auto tanh(SIMD<T, Size> x) noexcept {
        return SIMD<T, Size>(Physica::tanh(x.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> lncosh(SIMD<T, Size> x) noexcept {
        using Pack = SIMD<T, Size>;
        const auto x1 = abs(x);
        return x1 + ln1p(exp(-x1 * Pack(2))) - Pack(M_LN2);
    }

    template<Scalar T, int Size>
    auto round(SIMD<T, Size> x) noexcept -> SIMD<T, Size> {
        return Physica::round(x.toMachine());
    }
}

#include "Math/Exp.h"
