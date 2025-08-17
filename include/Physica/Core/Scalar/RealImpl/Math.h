/*
 * Copyright 2019-2025 Weibo He.
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

#include "Float64.h"

namespace Physica {
    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> unit(const Real<Prec>& x) noexcept {
        return Real<Prec>(x.isNegative() ? -1 : 1);
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> abs(const Real<Prec>& x) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::fabsf(x.toMachine()));
        else
            return Real<Prec>(::fabs(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> relu(const Real<Prec>& x) noexcept {
        return x.isPositive() ? x : Real<Prec>(0);
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> square(const Real<Prec>& x) noexcept {
        return x * x;
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> reciprocal(const Real<Prec>& x) noexcept {
        assert(!x.isZero() && "[Error]: Divide by zero");
        return Real<Prec>(1) / x;
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> sqrt(const Real<Prec>& x) noexcept {
        assert(!x.isNegative());
        return Real<Prec>(std::sqrt(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> cbrt(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::cbrt(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> ln(const Real<Prec>& x) noexcept {
        assert(x.isPositive() && "[Error]: Invalid param");
        return Real<Prec>(std::log(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> ln1p(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::log1p(x.toMachine()));
    }
    /**
     * \return log_a n
     */
    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> log(const Real<Prec>& x, const Real<Prec>& a) noexcept {
        return ln(x) / ln(a);
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> exp(const Real<Prec>& x) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::expf(x.toMachine()));
        else {
            static_assert(Prec == Float64, "[Error]: Unexpected type");
            return Real<Prec>(::exp(x.toMachine()));
        }
    }

    template<FloatPrec Prec>
    Real<Prec> pow(const Real<Prec>& x, const Real<Prec>& n) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::powf(x.toMachine(), n.toMachine()));
        else {
            static_assert(Prec == Float64, "[Error]: Unexpected type");
            return Real<Prec>(::pow(x.toMachine(), n.toMachine()));
        }
    }
    /*!
     * Ignoring error. If x is a float number, use floor() first.
     * \param x must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<FloatPrec Prec>
    Real<Prec> factorial(const Real<Prec>& x) noexcept {
        using FloatType = decltype(x.toMachine());
        const auto trivial = x.toMachine();
        FloatType temp = 1;
        FloatType result = 0;
        while(temp < trivial) {
            temp += 1;
            result *= temp;
        }
        return Real<Prec>(result);
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> cos(const Real<Prec>& x) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::cosf(x.toMachine()));
        else
            return Real<Prec>(::cos(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> cospi(const Real<Prec>& x) noexcept {
    #ifdef __CUDA_ARCH__
        if constexpr (Prec == Float32)
            return Real<Prec>(::cospif(x.toMachine()));
        else
            return Real<Prec>(::cospi(x.toMachine()));
    #else
        return cos(MathConst<Real<Prec>>::pi * x);
    #endif
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> sin(const Real<Prec>& x) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::sinf(x.toMachine()));
        else
            return Real<Prec>(::sin(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> sinpi(const Real<Prec>& x) noexcept {
    #ifdef __CUDA_ARCH__
        if constexpr (Prec == Float32)
            return Real<Prec>(::sinpif(x.toMachine()));
        else
            return Real<Prec>(::sinpi(x.toMachine()));
    #else
        return sin(MathConst<Real<Prec>>::pi * x);
    #endif
    }

    template<FloatPrec Prec>
    __host__ __device__ void sincos(Real<Prec> x, Real<Prec>& __restrict sin_result, Real<Prec>& __restrict cos_result) noexcept {
        if constexpr (Prec == Float32)
            ::sincosf(x.toMachine(), (float*)&sin_result, (float*)&cos_result);
        else
            ::sincos(x.toMachine(), (double*)&sin_result, (double*)&cos_result);
    }

    template<FloatPrec Prec>
    __host__ __device__ void sincospi(Real<Prec> x, Real<Prec>& __restrict sin_result, Real<Prec>& __restrict cos_result) noexcept {
    #ifdef __CUDA_ARCH__
        if constexpr (Prec == Float32)
            ::sincospif(x.toMachine(), (float*)&sin_result, (float*)&cos_result);
        else
            ::sincospi(x.toMachine(), (double*)&sin_result, (double*)&cos_result);
    #else
        return sincos(MathConst<Real<Prec>>::pi * x, sin_result, cos_result);
    #endif
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> tan(const Real<Prec>& x) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::tanf(x.toMachine()));
        else
            return Real<Prec>(::tan(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> tanpi(const Real<Prec>& x) noexcept {
        Real<Prec> sinpi_, cospi_;
        sincospi(x, sinpi_, cospi_);
        return sinpi_ / cospi_;
    }

    template<FloatPrec Prec>
    Real<Prec> sec(const Real<Prec>& x) noexcept {
        return reciprocal(cos(x));
    }

    template<FloatPrec Prec>
    Real<Prec> secpi(const Real<Prec>& x) noexcept {
        return reciprocal(cospi(x));
    }

    template<FloatPrec Prec>
    Real<Prec> csc(const Real<Prec>& x) noexcept {
        return reciprocal(sin(x));
    }

    template<FloatPrec Prec>
    Real<Prec> cscpi(const Real<Prec>& x) noexcept {
        return reciprocal(sinpi(x));
    }

    template<FloatPrec Prec>
    Real<Prec> cot(const Real<Prec>& x) noexcept {
        Real<Prec> sin_, cos_;
        sincos(x, sin_, cos_);
        return cos_ / sin_;
    }

    template<FloatPrec Prec>
    Real<Prec> cotpi(const Real<Prec>& x) noexcept {
        Real<Prec> sinpi_, cospi_;
        sincospi(x, sinpi_, cospi_);
        return cospi_ / sinpi_;
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> arccos(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::acos(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> arcsin(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::asin(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> arctan(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::atan(x.toMachine()));
    }

    template<FloatPrec Prec>
    Real<Prec> arcsec(const Real<Prec>& x) noexcept {
        return arccos(reciprocal(x));
    }

    template<FloatPrec Prec>
    Real<Prec> arccsc(const Real<Prec>& x) noexcept {
        return arcsin(reciprocal(x));
    }

    template<FloatPrec Prec>
    Real<Prec> arccot(const Real<Prec>& x) noexcept {
        return arctan(reciprocal(x));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> cosh(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::cosh(x.toMachine()));
    }


    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> sinh(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::sinh(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> tanh(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::tanh(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> sech(const Real<Prec>& x) noexcept {
        return Real<Prec>(1 / std::cosh(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> csch(const Real<Prec>& x) noexcept {
        return Real<Prec>(1 / std::sinh(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> coth(const Real<Prec>& x) noexcept {
        return Real<Prec>(1 / std::tanh(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> arccosh(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::acosh(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> arcsinh(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::asinh(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> arctanh(const Real<Prec>& x) noexcept {
        return Real<Prec>(std::atanh(x.toMachine()));
    }

    template<FloatPrec Prec>
    Real<Prec> arcsech(const Real<Prec>& x) noexcept {
        return arccosh(reciprocal(x));
    }

    template<FloatPrec Prec>
    Real<Prec> arccsch(const Real<Prec>& x) noexcept {
        return arcsinh(reciprocal(x));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> arccoth(const Real<Prec>& x) noexcept {
        auto trivial = x.toMachine();
        return Real<Prec>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }

    template<FloatPrec Prec>
    Real<Prec> ln1pexp(const Real<Prec>& x) noexcept {
        return relu(x) + ln1p(exp(-abs(x)));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> lncosh(const Real<Prec>& x) noexcept {
        using T = Real<Prec>;
        const auto x1 = abs(x);
        return x1 + ln1p(exp(-x1 * T(2))) - T(M_LN2);
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> floor(const Real<Prec>& x) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::floorf(x.toMachine()));
        else
            return Real<Prec>(::floor(x.toMachine()));
    }

    template<FloatPrec Prec>
    __host__ __device__ Real<Prec> ceil(const Real<Prec>& x) noexcept {
        if constexpr (Prec == Float32)
            return Real<Prec>(::ceilf(x.toMachine()));
        else
            return Real<Prec>(::ceil(x.toMachine()));
    }
}
