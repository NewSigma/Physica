/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Scalar/Real.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
    #include <vectorclass/vectormath_exp.h>
    #include <vectorclass/vectormath_trig.h>
    #include <vectorclass/vectormath_hyp.h>
#pragma GCC diagnostic pop

using namespace Physica;

namespace Physica {
    template<FloatPrec Prec, int Size>
    auto ln(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(log(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto ln1p(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(log1p(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto pow(SIMD<Real<Prec>, Size> x, SIMD<Real<Prec>, Size> y) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(pow(x.toMachine(), y.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto sin(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        Internal::checkTrigonometricParam<false>(abs(x).max());
        return SIMD<Real<Prec>, Size>(VCL::sin(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto cos(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        Internal::checkTrigonometricParam<false>(abs(x).max());
        return SIMD<Real<Prec>, Size>(VCL::cos(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto sincos(SIMD<Real<Prec>, Size> x) noexcept -> std::pair<SIMD<Real<Prec>, Size>, SIMD<Real<Prec>, Size>> {
        Internal::checkTrigonometricParam<false>(abs(x).max());
        SIMD<Real<Prec>, Size> s, c;
        VCL::sincos(x.toMachine(), s.toMachine(), c.toMachine());
        return std::make_pair(s, c);
    }

    template<FloatPrec Prec, int Size>
    auto tan(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        Internal::checkTrigonometricParam<false>(abs(x).max());
        return SIMD<Real<Prec>, Size>(VCL::tan(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto arctan2(SIMD<Real<Prec>, Size> y, SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return VCL::atan2(y.toMachine(), x.toMachine());
    }

    template<FloatPrec Prec, int Size>
    auto cosh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(VCL::cosh(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto sinh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(VCL::sinh(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto tanh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(VCL::tanh(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto arcsinh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(VCL::asinh(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto arctanh(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(VCL::atanh(x.toMachine()));
    }

    template<FloatPrec Prec, int Size>
    auto cbrt(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        return SIMD<Real<Prec>, Size>(cbrt(x.toMachine()));
    }
}

#define PHYSICA_INSTANTIATE(Prec, Size) \
namespace Physica { \
    template auto ln(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto ln1p(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto pow(SIMD<Real<Prec>, Size>, SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto sin(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto cos(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto sincos(SIMD<Real<Prec>, Size>) noexcept -> std::pair<SIMD<Real<Prec>, Size>, SIMD<Real<Prec>, (Size)>>; \
    template auto tan(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto arctan2(SIMD<Real<Prec>, Size>, SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto cosh(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto sinh(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto tanh(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto arcsinh(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto arctanh(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
    template auto cbrt(SIMD<Real<Prec>, Size>) noexcept -> SIMD<Real<Prec>, Size>; \
}

PHYSICA_INSTANTIATE(Float32, 4)
PHYSICA_INSTANTIATE(Float32, 8)
PHYSICA_INSTANTIATE(Float32, 16)
PHYSICA_INSTANTIATE(Float64, 2)
PHYSICA_INSTANTIATE(Float64, 4)
PHYSICA_INSTANTIATE(Float64, 8)
