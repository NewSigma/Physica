/*
 * Copyright 2020-2026 Weibo He.
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

#include "../Complex.h"

namespace Physica {
    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> unit(const Complex<T>& x) noexcept {
        const T temp = x.norm();
        if (temp.isZero())
            return T(1);
        return x * reciprocal(temp);
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> fma(Complex<T> x, Complex<T> y, Complex<T> z) noexcept {
        T re = fma(x.imag(), -y.imag(), fma(x.real(), y.real(), z.real()));
        T im = fma(x.imag(), y.real(), fma(x.real(), y.imag(), z.imag()));
        return Complex<T>(re, im);
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> fma(Complex<T> x, T y, Complex<T> z) noexcept {
        return Complex<T>(fma(x.real(), y, z.real()), fma(x.imag(), y, z.imag()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ T abs(const Complex<T>& c) noexcept {
        return c.norm();
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> square(const Complex<T>& c) noexcept {
        const auto& re = c.real();
        const auto& im = c.imag();
        return Complex<T>(fma(re, re, -square(im)), (re * im) * T(2));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> reciprocal(const Complex<T>& c) noexcept {
        assert(!c.isSubNormal() && "[Error]: Division overflow");
        const auto& real = c.real();
        const auto& imag = c.imag();
        const auto divisor = reciprocal(square(real) + square(imag));
        return Complex<T>(real * divisor, -imag * divisor);
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. Numerical Recipes(3rd edition)[M]. London: Cambridge University Press, 2007:226
     */
    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> sqrt(const Complex<T>& c) noexcept {
        using ResultType = Complex<T>;
        if (c.isZero())
            return ResultType(0);
        using Tr = T::RealType;
        const Tr abs_real = abs(c.real());
        const Tr w = sqrt((abs_real + c.norm()) * T(0.5));
        const Tr v = c.imag() / w * T(0.5);
        if (!c.real().isNegative())
            return ResultType(w, v);
        return ResultType(abs(v), c.imag().isNegative() ? -w : w);
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> ln(const Complex<T>& c) noexcept {
        assert(!c.isZero());
        if constexpr (IsHost())
            return Complex<T>(std::log(c.toMachine()));
        else {
        #ifdef PHYSICA_CUDA
            return Complex<T>(thrust::log(c.toThrust()));
        #else
            unreachable();
        #endif
        }
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> ln1p(const Complex<T>& c) noexcept {
        return ln(T(1) + c);
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> exp(const Complex<T>& c) noexcept {
        return Complex<T>(std::exp(c.toMachine()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> expm1(const Complex<T>& c) noexcept {
        return exp(c) - T(1);
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> cos(const Complex<T>& c) noexcept {
        return Complex<T>(std::cos(c.toMachine()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> sin(const Complex<T>& c) noexcept {
        return Complex<T>(std::sin(c.toMachine()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> tan(const Complex<T>& c) noexcept {
        return Complex<T>(std::tan(c.toMachine()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> sec(const Complex<T>& c) noexcept {
        return reciprocal(cos(c));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> csc(const Complex<T>& c) noexcept {
        return reciprocal(sin(c));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> cot(const Complex<T>& c) noexcept {
        return reciprocal(tan(c));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> cosh(const Complex<T>& c) noexcept {
        return Complex<T>(std::cosh(c.toMachine()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> sinh(const Complex<T>& c) noexcept {
        return Complex<T>(std::sinh(c.toMachine()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> tanh(const Complex<T>& c) noexcept {
        return Complex<T>(std::tanh(c.toMachine()));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> sech(const Complex<T>& c) noexcept {
        return reciprocal(cosh(c));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> csch(const Complex<T>& c) noexcept {
        return reciprocal(sinh(c));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> coth(const Complex<T>& c) noexcept {
        return reciprocal(tanh(c));
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> lncosh(const Complex<T>& c) noexcept {
        const T abs_real = abs(c.real());
        const T norm1 = exp(T(-2) * abs_real);
        auto [sine, cosine] = sincos(c.real().isPositive() ? c.imag() : -c.imag());
        const auto temp = Complex<T>(fma(norm1, cosine, cosine), fma(-norm1, sine, sine)) * T(0.5);
        return abs_real + ln(temp);
    }

    template<Scalar T>
    [[nodiscard]] __host__ __device__ Complex<T> softplus(const Complex<T>& c) noexcept {
        if (c.real().isPositive()) {
            T norm = c.norm();
            return norm + ln(exp(-norm) + exp(c - norm));
        }
        return ln1p(exp(c));
    }
}
