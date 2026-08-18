/*
 * Copyright 2025 Weibo He.
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

#include "../Math.h"

namespace Physica {
    /**
     * Adapted from [1], assuming no-signed-zero and finite-math
     * This function does not produce denormals $^{[1]}$.
     *
     * \tparam
     * Base: 0 for exp, 1 for 0.5 * exp, 2 for pow(2, x), 10 for pow(10, x)
     *
     * Reference:
     * [1] vectorclass2; https://github.com/vectorclass/version2/blob/master/vectormath_exp.h
     */
    namespace Internal {

        template<Packet Pack, bool ExpM1, int Base>
        Pack exp_float32(const Pack initial_x) {
            // Taylor coefficients
            // Not using minimax approximation because we prioritize precision close to x = 0
            constexpr float P0 = 1.0 / 2;
            constexpr float P1 = 1.0 / 6;
            constexpr float P2 = 1.0 / 24;
            constexpr float P3 = 1.0 / 120;
            constexpr float P4 = 1.0 / 720;
            constexpr float P5 = 1.0 / 5040;
            // The lower limit of x is slightly more restrictive than the upper limit.
            // We are specifying the lower limit, except for Base = 1 because it is only for positive x in hyperbolic functions
            constexpr float32 max_x = []() consteval static noexcept -> float32 {
                switch (Base) {
                case 0:
                    return 87.33655F;
                case 1:
                    return 89.415985F;
                case 2:
                    return 126.0F;
                case 10:
                    return 37.92978F;
                default:
                    unreachable();
                };
            }();

            Pack x;
            Pack r;
            if constexpr (Base <= 1) {
                constexpr float ln2_hi = 0.693359375F;
                constexpr float ln2_lo = -2.12194440E-4F;
                r = round(initial_x * std::numbers::log2e_v<float>);
                x = initial_x;
                x = nmul_add(r, Pack(ln2_hi), x); // subtraction in two steps for higher precision
                x = nmul_add(r, Pack(ln2_lo), x);
                if constexpr (Base == 1)
                    r--; // 0.5 * exp(x)
            }
            else if constexpr (Base == 2) {
                r = round(initial_x);
                x = (initial_x - r) * std::numbers::ln2_v<float>;
            }
            else {
                constexpr float log10_2_hi = 0.301025391F;
                constexpr float log10_2_lo = 4.60503907E-6F;
                r = round(initial_x * (std::numbers::log2e_v<float> * std::numbers::ln10_v<float>));
                x = initial_x;
                x = nmul_add(r, Pack(log10_2_hi), x);
                x = nmul_add(r, Pack(log10_2_lo), x);
                x = x * std::numbers::ln10_v<float>;
            }

            Pack z = polynomial_5(x.toMachine(), P0, P1, P2, P3, P4, P5);
            const Pack n2 = vm_pow2n(r.toMachine());
            z = fma(z, x * x, x);
            z = ExpM1 ? fma(z, n2, n2 - Pack(1)) : ((z + Pack(1)) * n2);

            const auto inRange  = abs(initial_x) < Pack(max_x);
            const bool overflow = !inRange.horizontal_and();
            if (overflow) [[unlikely]] {
                const Pack z0 = Pack::select(initial_x.isPositive(), Pack::inf(), Pack(ExpM1 ? -1 : 0));
                z = Pack::select(inRange, z, z0);
            }
            return z;
        }

        template<Packet Pack, bool ExpM1, int Base>
        Pack exp_float64(const Pack initial_x) noexcept {
            constexpr double P2  = 1.0 / 2;
            constexpr double P3  = 1.0 / 6;
            constexpr double P4  = 1.0 / 24;
            constexpr double P5  = 1.0 / 120;
            constexpr double P6  = 1.0 / 720;
            constexpr double P7  = 1.0 / 5040;
            constexpr double P8  = 1.0 / 40320;
            constexpr double P9  = 1.0 / 362880;
            constexpr double P10 = 1.0 / 3628800;
            constexpr double P11 = 1.0 / 39916800;
            constexpr double P12 = 1.0 / 479001600;
            constexpr double P13 = 1.0 / 6227020800;
            constexpr float64 max_x = []() consteval static noexcept -> float64 {
                switch (Base) {
                case 0:
                    return 708.3964185322641;
                case 1:
                    return 709.782712893384;
                case 2:
                    return 1022.0;
                case 10:
                    return 307.6526555685888;
                default:
                    unreachable();
                };
            }();

            Pack x;
            Pack r;
            if constexpr (Base <= 1) {
                constexpr double ln2_hi = 0.693145751953125;
                constexpr double ln2_lo = 1.42860682030941723212E-6;
                r = round(initial_x * std::numbers::log2e);
                x = initial_x;
                x = nmul_add(r, Pack(ln2_hi), x);
                x = nmul_add(r, Pack(ln2_lo), x);
                if constexpr (Base == 1)
                    r--;
            }
            else if constexpr (Base == 2) {
                r = round(initial_x);
                x = (initial_x - r) * std::numbers::ln2;
            }
            else {
                constexpr double log10_2_hi = 0.30102999554947019;
                constexpr double log10_2_lo = 1.1451100899212592E-10;
                r = round(initial_x * (std::numbers::log2e * std::numbers::ln10));
                x = initial_x;
                x = nmul_add(r, Pack(log10_2_hi), x);
                x = nmul_add(r, Pack(log10_2_lo), x);
                x *= std::numbers::ln10;
            }

            Pack z = polynomial_13m(x.toMachine(), P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13);
            const Pack n2 = vm_pow2n(r.toMachine());
            z = ExpM1 ? fma(z, n2, n2 - Pack(1)) : (z + Pack(1)) * n2;

            const auto inRange  = abs(initial_x) < Pack(max_x);
            const bool overflow = !inRange.horizontal_and();
            if (overflow) [[unlikely]] {
                const Pack z0 = Pack::select(initial_x.isPositive(), Pack::inf(), Pack(ExpM1 ? -1 : 0));
                z = Pack::select(inRange, z, z0);
            }
            return z;
        }
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto exp(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        if constexpr (Prec == Float32)
            return Internal::exp_float32<SIMD<Real<Prec>, Size>, false, 0>(x);
        else
            return Internal::exp_float64<SIMD<Real<Prec>, Size>, false, 0>(x);
    }

    template<FloatPrec Prec, int Size>
    [[nodiscard]] auto expm1(SIMD<Real<Prec>, Size> x) noexcept -> SIMD<Real<Prec>, Size> {
        if constexpr (Prec == Float32)
            return Internal::exp_float32<SIMD<Real<Prec>, Size>, true, 0>(x);
        else
            return Internal::exp_float64<SIMD<Real<Prec>, Size>, true, 0>(x);
    }
}
