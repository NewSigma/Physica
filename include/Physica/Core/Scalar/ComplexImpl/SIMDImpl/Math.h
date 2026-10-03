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

#include "Physica/Core/Scalar/ComplexImpl/SIMD.h"

namespace Physica {
    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> unit(const SIMD<Complex<T>, Size> x) noexcept {
        Array<Complex<T>, Size> buffer;
        for (int i = 0; i < Size; ++i)
            buffer[i] = unit(x[i]);

        SIMD<Complex<T>, Size> result;
        result.load(buffer.data());
        return result;
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> fma(const SIMD<Complex<T>, Size> a, const SIMD<Complex<T>, Size> b, const SIMD<Complex<T>, Size> c) noexcept {
        return a * b + c;
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> fma(const SIMD<Complex<T>, Size> a, const SIMD<T, Size> b, const SIMD<Complex<T>, Size> c) noexcept {
        using FullRealType = SIMD<Complex<T>, Size>::FullRealType;
        return SIMD<Complex<T>, Size>::asComplex(fma(a.asReal(), FullRealType(b, b).gatherRealImag(), c.asReal()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> fma(const SIMD<T, Size> a, const SIMD<Complex<T>, Size> b, const SIMD<Complex<T>, Size> c) noexcept {
        return fma(b, a, c);
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size * 2> abs(const SIMD<Complex<T>, Size> x) noexcept {
        return sqrt(x.squaredNorm());
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> square(const SIMD<Complex<T>, Size> x) noexcept {
        return x * x;
    }
    /**
     * Mirroring the scalar path
     *
     * References:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. Numerical Recipes(3rd edition)[M]. London: Cambridge University Press, 2007:226
     */
    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> sqrt(const SIMD<Complex<T>, Size> x) noexcept {
        using ResultType = SIMD<Complex<T>, Size>;
        using RealPack = SIMD<T, Size * 2>;
        const auto [real, imag] = x.makeFullRealImag();

        const RealPack abs_real = abs(real);
        const RealPack norm = sqrt(fma(real, real, square(imag)));
        const RealPack w = sqrt((abs_real + norm) * T(0.5));
        const RealPack v = imag / w * T(0.5);

        const auto isNeg = real.isNegative();
        const RealPack result_real = RealPack::select(isNeg, abs(v), w);
        const RealPack result_imag = RealPack::select(isNeg, RealPack::select(imag.isNegative(), -w, w), v);

        RealPack result;
        if constexpr (Size == 1)
            result = RealPack::template blend<0, 2>(result_real, result_imag);
        else if constexpr (Size == 2)
            result = RealPack::template blend<0, 4, 2, 6>(result_real, result_imag);
        else if constexpr (Size == 4)
            result = RealPack::template blend<0, 8, 2, 10, 4, 12, 6, 14>(result_real, result_imag);
        else {
            static_assert(Size == 8, "[Error]: Unexpected size");
            result = RealPack::template blend<0, 16, 2, 18, 4, 20, 6, 22, 8, 24, 10, 26, 12, 28, 14, 30>(result_real, result_imag);
        }
        return ResultType::asComplex(RealPack::select(norm.isZero(), RealPack(0), result));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> exp(const SIMD<Complex<T>, Size> x) noexcept {
        using ResultType = SIMD<Complex<T>, Size>;
        using RealType = ResultType::RealType;
        if constexpr (ResultType::isSeparatable) {
            using FullRealType = ResultType::FullRealType;
            const RealType factor = exp(x.real());
            auto [s, c] = sincos(x.imag());
            s *= factor;
            c *= factor;
            return ResultType::asComplex(FullRealType(c, s).scatterRealImag());
        }
        else {
            const auto [re, im] = x.makeFullRealImag();
            const RealType factor = exp(re);
            const auto [s, c] = sincos(im);
            RealType cs;
            if constexpr (T::Prec == Float32)
                cs = RealType::template blend<0, 4, 3, 7>(c, s);
            else
                cs = RealType::template blend<0, 2>(c, s);
            return ResultType::asComplex(factor * cs);
        }
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> ln(const SIMD<Complex<T>, Size> x) noexcept {
        using ResultType = SIMD<Complex<T>, Size>;
        using RealType = ResultType::RealType;
        if constexpr (ResultType::isSeparatable) {
            using FullRealType = ResultType::FullRealType;
            const auto x1 = ResultType::asComplex(x.asReal());
            auto x_re = x1.real();
            auto x_im = x1.imag();

            const auto factor = reciprocal(std::max(abs(x_re), abs(x_im)));
            const auto lnF = ln(factor);
            x_re *= factor;
            x_im *= factor;

            const auto re = mul_sub(ln(fma(x_re, x_re, square(x_im))), RealType(0.5), RealType(lnF));
            const auto im = arctan2(x_im, x_re);
            return ResultType::asComplex(FullRealType(re, im).scatterRealImag());
        }
        else {
            const auto x1 = abs(x.asReal());
            const auto factor = reciprocal(std::max(x1, x1.swapRealImag()));
            const auto lnF = ln(factor);

            const auto x2 = ResultType::asComplex(x.asReal() * factor);
            const auto re = mul_sub(ln(x2.squaredNorm()), RealType(0.5), RealType(lnF));
            const auto [x2re, x2im] = x2.makeFullRealImag();
            const auto im = arctan2(x2im, x2re);

            RealType result;
            if constexpr (T::Prec == Float32)
                result = RealType::template blend<0, 4, 3, 7>(re, im);
            else
                result = RealType::template blend<0, 2>(re, im);
            return ResultType::asComplex(result);
        }
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> ln1p(const SIMD<Complex<T>, Size> x) noexcept {
        return ln(T(1) + x);
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> expm1(const SIMD<Complex<T>, Size> x) noexcept {
        return exp(x) - SIMD<Complex<T>, Size>(T(1));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> cosh(const SIMD<Complex<T>, Size> x) noexcept {
        return (exp(x) + exp(-x)) * T(0.5);
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> sinh(const SIMD<Complex<T>, Size> x) noexcept {
        return (exp(x) - exp(-x)) * T(0.5);
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> tanh(const SIMD<Complex<T>, Size> x) noexcept {
        using ResultType = SIMD<Complex<T>, Size>;
        std::array<Complex<T>, Size> arr;
        for (int i = 0; i < Size; ++i)
            arr[i] = tanh(x[i]);
        ResultType result{};
        result.load(arr.data());
        return result;
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> sech(const SIMD<Complex<T>, Size> x) noexcept {
        return T(0.5) / (exp(x) + exp(-x));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> lncosh(const SIMD<Complex<T>, Size> x) noexcept {
        using ResultType = SIMD<Complex<T>, Size>;
        using RealType = ResultType::RealType;
        if constexpr (ResultType::isSeparatable) {
            using FullRealType = ResultType::FullRealType;
            const RealType re = x.real();
            const RealType im = x.imag();
            const RealType abs_real = abs(x.real());
            const RealType norm1 = exp(T(-2) * abs_real);
            const auto [s, c] = sincos(RealType::select(re.isPositive(), im, -im));
            const auto temp = ResultType::asComplex(FullRealType(fma(norm1, c, c), nmul_add(norm1, s, s)).scatterRealImag() * T(0.5));
            return ResultType::asComplex(FullRealType(abs_real, RealType(0)).scatterRealImag()) + ln(temp);
        }
        else {
            const auto [re, im] = x.makeFullRealImag();
            const RealType abs_real = abs(re);
            const RealType norm1 = exp(T(-2) * abs_real);
            const auto [s, c] = sincos(RealType::select(re.isPositive(), -im, im));
            RealType cs;
            if constexpr (T::Prec == Float32)
                cs = RealType::template blend<0, 4, 3, 7>(c, s);
            else
                cs = RealType::template blend<0, 2>(c, s);
            const auto temp = ResultType::asComplex(mul_addsub(norm1, cs, -cs) * T(0.5));
            return ResultType::asComplex(ResultType::asComplex(abs_real).real()) + ln(temp);
        }
    }
}
