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
    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> abs(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::fabsf(x.toMachine()));
        else
            return Real<Option>(::fabs(x.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> relu(const Real<Option>& x) noexcept {
        return x.isPositive() ? x : Real<Option>(0);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> square(const Real<Option>& x) noexcept {
        return x * x;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> reciprocal(const Real<Option>& x) noexcept {
        if constexpr (Option != FloatMP)
            assert(!x.isZero() && "[Error]: Divide by zero");
        return Real<Option>(1) / x;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> sqrt(const Real<Option>& x) noexcept {
        assert(!x.isNegative());
        if constexpr (Option == Float32)
            return Real<Option>(::sqrtf(x.toMachine()));
        else {
            static_assert(Option == Float64, "[Error]: Unexpected type");
            return Real<Option>(::sqrt(x.toMachine()));
        }
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> cbrt(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::cbrtf(x.toMachine()));
        else {
            static_assert(Option == Float64, "[Error]: Unexpected type");
            return Real<Option>(::cbrt(x.toMachine()));
        }
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ln(const Real<Option>& x) noexcept {
        assert(x.isPositive() && "[Error]: Invalid param");
        return Real<Option>(std::log(x.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ln1p(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::log1pf(x.toMachine()));
        else {
            static_assert(Option == Float64, "[Error]: Unexpected type");
            return Real<Option>(::log1p(x.toMachine()));
        }
    }
    /**
     * \return log_a n
     */
    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> log(const Real<Option>& x, const Real<Option>& a) noexcept {
        return ln(x) / ln(a);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> exp(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::expf(x.toMachine()));
        else {
            static_assert(Option == Float64, "[Error]: Unexpected type");
            return Real<Option>(::exp(x.toMachine()));
        }
    }

    template<ScalarOption Option>
    Real<Option> pow(const Real<Option>& x, const Real<Option>& n) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::powf(x.toMachine(), n.toMachine()));
        else {
            static_assert(Option == Float64, "[Error]: Unexpected type");
            return Real<Option>(::pow(x.toMachine(), n.toMachine()));
        }
    }
    /*!
     * Ignoring error. If x is a float number, use floor() first.
     * \param x must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<ScalarOption Option>
    Real<Option> factorial(const Real<Option>& x) noexcept {
        using FloatType = decltype(x.toMachine());
        const auto trivial = x.toMachine();
        FloatType temp = 1;
        FloatType result = 0;
        while(temp < trivial) {
            temp += 1;
            result *= temp;
        }
        return Real<Option>(result);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> cos(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::cosf(x.toMachine()));
        else
            return Real<Option>(::cos(x.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> sin(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::sinf(x.toMachine()));
        else
            return Real<Option>(::sin(x.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline void sincos(Real<Option> x, Real<Option>& sin_result, Real<Option>& cos_result) noexcept {
        using MachineType = Real<Option>::MachineType;
        MachineType sin_temp, cos_temp;
        if constexpr (Option == Double)
            ::sincos(x.toMachine(), (double*)&sin_temp, (double*)&cos_temp);
        else
            ::sincosf(x.toMachine(), (float*)&sin_temp, (float*)&cos_temp);
        sin_result = sin_temp;
        cos_result = cos_temp;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> tan(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::tanf(x.toMachine()));
        else
            return Real<Option>(::tan(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> sec(const Real<Option>& x) noexcept {
        return reciprocal(cos(x));
    }

    template<ScalarOption Option>
    Real<Option> csc(const Real<Option>& x) noexcept {
        return reciprocal(sin(x));
    }

    template<ScalarOption Option>
    Real<Option> cot(const Real<Option>& x) noexcept {
        return cos(x) / sin(x);
    }

    template<ScalarOption Option>
    Real<Option> arccos(const Real<Option>& x) noexcept {
        return Real<Option>(std::acos(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsin(const Real<Option>& x) noexcept {
        return Real<Option>(std::asin(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arctan(const Real<Option>& x) noexcept {
        return Real<Option>(std::atan(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsec(const Real<Option>& x) noexcept {
        return arccos(reciprocal(x));
    }

    template<ScalarOption Option>
    Real<Option> arccsc(const Real<Option>& x) noexcept {
        return arcsin(reciprocal(x));
    }

    template<ScalarOption Option>
    Real<Option> arccot(const Real<Option>& x) noexcept {
        return arctan(reciprocal(x));
    }

    template<ScalarOption Option>
    Real<Option> cosh(const Real<Option>& x) noexcept {
        return Real<Option>(std::cosh(x.toMachine()));
    }


    template<ScalarOption Option>
    Real<Option> sinh(const Real<Option>& x) noexcept {
        return Real<Option>(std::sinh(x.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ Real<Option> tanh(const Real<Option>& x) noexcept {
        return Real<Option>(std::tanh(x.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ Real<Option> sech(const Real<Option>& x) noexcept {
        return Real<Option>(1 / std::cosh(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> csch(const Real<Option>& x) noexcept {
        return Real<Option>(1 / std::sinh(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> coth(const Real<Option>& x) noexcept {
        return Real<Option>(1 / std::tanh(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arccosh(const Real<Option>& x) noexcept {
        return Real<Option>(std::acosh(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsinh(const Real<Option>& x) noexcept {
        return Real<Option>(std::asinh(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arctanh(const Real<Option>& x) noexcept {
        return Real<Option>(std::atanh(x.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsech(const Real<Option>& x) noexcept {
        return arccosh(reciprocal(x));
    }

    template<ScalarOption Option>
    Real<Option> arccsch(const Real<Option>& x) noexcept {
        return arcsinh(reciprocal(x));
    }

    template<ScalarOption Option>
    Real<Option> arccoth(const Real<Option>& x) noexcept {
        auto trivial = x.toMachine();
        return Real<Option>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }

    template<ScalarOption Option>
    Real<Option> ln1pexp(const Real<Option>& x) noexcept {
        return relu(x) + ln1p(exp(-abs(x)));
    }

    template<ScalarOption Option>
    __host__ __device__ Real<Option> lncosh(const Real<Option>& x) noexcept {
        using T = Real<Option>;
        const auto x1 = abs(x);
        return x1 + ln1p(exp(-x1 * T(2))) - T(M_LN2);
    }

    template<ScalarOption Option>
    __host__ __device__ Real<Option> floor(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::floorf(x.toMachine()));
        else
            return Real<Option>(::floor(x.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ceil(const Real<Option>& x) noexcept {
        if constexpr (Option == Float32)
            return Real<Option>(::ceilf(x.toMachine()));
        else
            return Real<Option>(::ceil(x.toMachine()));
    }
}
