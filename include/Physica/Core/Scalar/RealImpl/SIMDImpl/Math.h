/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

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
    [[nodiscard]] PHYSICA_API auto cbrt(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto ln(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto ln1p(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto exp(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto expm1(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto pow(SIMD<Real<Prec>, Size> x, SIMD<Real<Prec>, Size> y) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto sin(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto cos(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto sincos(SIMD<Real<Prec>, Size> x) noexcept -> std::pair<SIMD<Real<Prec>, Size>, SIMD<Real<Prec>, Size>>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto tan(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto sec(SIMD<Real<Prec>, Size> x) noexcept {
        return reciprocal(cos(x));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto arctan2(SIMD<Real<Prec>, Size> y, SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto cosh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto sinh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto tanh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto sech(SIMD<Real<Prec>, Size> x) noexcept {
        return reciprocal(cosh(x));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto arcsinh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] PHYSICA_API auto arctanh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size>;

    template<FloatPrec Prec, int Size>
    [[nodiscard]] SIMD<Real<Prec>, Size> lncosh(SIMD<Real<Prec>, Size> x) noexcept {
        using Pack = SIMD<Real<Prec>, Size>;
        const auto x1 = abs(x);
        const auto Overflow = x1 > Pack(Internal::lncosh_threshold<Prec>());
        if (!Overflow.horizontal_or())
            return Pack(0.5) * ln1p(square(sinh(x1)));

        const Pack large = x1 + ln1p(exp(-x1 * Pack(2))) - Pack(M_LN2);
        if (Overflow.horizontal_and())
            return large;

        const Pack small = Pack(0.5) * ln1p(square(sinh(x1)));
        return Pack::select(Overflow, large, small);
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto sigmoid(SIMD<Real<Prec>, Size> x) noexcept {
        using Pack = SIMD<Real<Prec>, Size>;
        const auto y = exp(-abs(x));
        return Pack::select(x.isPositive(), reciprocal(Pack(1) + y), y / (Pack(1) + y));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto softplus(SIMD<Real<Prec>, Size> x) noexcept {
        return relu(x) + ln1p(exp(-abs(x)));
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto round(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return VCL::round(x.toMachine());
    }
}
