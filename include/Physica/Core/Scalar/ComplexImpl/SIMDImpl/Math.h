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
        const auto [re, im] = a.makeFullRealImag();
        return SIMD<Complex<T>, Size>::asComplex(mul_addsub(re, b.asReal(), mul_addsub(im, b.swapRealImag(), c.asReal())));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> fma(const SIMD<Complex<T>, Size> a, const SIMD<T, Size> b, const SIMD<Complex<T>, Size> c) noexcept {
        using FullRealType = SIMD<Complex<T>, Size>::FullRealType;
        return SIMD<Complex<T>, Size>::asComplex(fma(a.asReal(), FullRealType(b, b).gatherRealImag(), c.asReal()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size * 2> abs(const SIMD<Complex<T>, Size> x) noexcept {
        return sqrt(x.squaredNorm());
    }
    /**
     * References:
     * [1] add-on; https://github.com/vectorclass/add-on
     */
    template<Scalar T, int Size>
    [[nodiscard]] SIMD<Complex<T>, Size> sqrt(const SIMD<Complex<T>, Size> x) noexcept {
        using RealPack = SIMD<T, Size * 2>;
        const RealPack x1 = x.asReal();
        const RealPack t1 = x1 * x1;
        RealPack t2;
        if constexpr (T::Prec == Float32)
            t2 = t1.template shuffle<1, 0, 3, 2>();
        else {
            static_assert(T::Prec == Float64, "[Error]: Not implemented");
            if constexpr (Size == 1)
                t2 = t1.template shuffle<1, 0>();
            else if constexpr (Size == 2)
                t2 = t1.template shuffle<1, 0, 1, 0>();
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                t2 = t1.template shuffle<1, 0, 1, 0, 1, 0, 1, 0>();
            }
        }

        const RealPack t3 = sqrt(t1 + t2);
        RealPack t4;
        if constexpr (T::Prec == Float32)
            t4 = x1.template shuffle<0, 0, 2, 2>();
        else {
            static_assert(T::Prec == Float64, "[Error]: Not implemented");
            if constexpr (Size == 1)
                t4 = x1.template shuffle<0, 0>();
            else if constexpr (Size == 2)
                t4 = x1.template shuffle<0, 0, 0, 0>();
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                t4 = x1.template shuffle<0, 0, 0, 0, 0, 0, 0, 0>();
            }
        }

        RealPack signbit;
        if constexpr (Size == 1)
            signbit = RealPack::template makeSignBits<0, 1>();
        else if constexpr (Size == 2)
            signbit = RealPack::template makeSignBits<0, 1, 0, 1>();
        else if constexpr (Size == 4)
            signbit = RealPack::template makeSignBits<0, 1, 0, 1, 0, 1, 0, 1>();
        else {
            static_assert(Size == 8, "[Error]: Unexpected size");
            signbit = RealPack::template makeSignBits<0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1>();
        }
        const RealPack t5 = t3 + (t4 ^ signbit);
        const RealPack t6 = sqrt(t5 * T(0.5));
        const RealPack result = t6 ^ (x1 & signbit);
        return SIMD<Complex<T>, Size>::asComplex(result);
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
    [[nodiscard]] SIMD<Complex<T>, Size> tanh(const SIMD<Complex<T>, Size> x) noexcept {
        using ResultType = SIMD<Complex<T>, Size>;
        std::array<Complex<T>, Size> arr;
        for (int i = 0; i < Size; ++i)
            arr[i] = tanh(x[i]);
        ResultType result{};
        result.load(arr);
        return result;
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
