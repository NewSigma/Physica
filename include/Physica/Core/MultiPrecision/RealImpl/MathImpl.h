/*
 * Copyright 2020-2024 Weibo He.
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

namespace Physica::Core {
    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> abs(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::fabsf(s.toMachine()));
        else
            return Real<Option>(::fabs(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> relu(const Real<Option>& s) noexcept {
        return s.isPositive() ? s : Real<Option>(0);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> square(const Real<Option>& s) noexcept {
        return s * s;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> reciprocal(const Real<Option>& s) noexcept {
        assert(!s.isZero() && "[Error]: Divide by zero");
        return Real<Option>(1) / s;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> sqrt(const Real<Option>& s) noexcept {
        assert(!s.isNegative());
        if constexpr (Option == Float)
            return Real<Option>(::sqrtf(s.toMachine()));
        else
            return Real<Option>(::sqrt(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> cbrt(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::cbrtf(s.toMachine()));
        else
            return Real<Option>(::cbrt(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ln(const Real<Option>& s) noexcept {
        assert(s.isPositive() && "[Error]: Invalid param");
        if constexpr (Option == Float)
            return Real<Option>(::logf(s.toMachine()));
        else
            return Real<Option>(::log(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ln1p(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::log1pf(s.toMachine()));
        else
            return Real<Option>(::log1p(s.toMachine()));
    }
    /**
     * \return log_a n
     */
    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> log(const Real<Option>& s, const Real<Option>& a) noexcept {
        return ln(s) / ln(a);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> exp(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::expf(s.toMachine()));
        else
            return Real<Option>(::exp(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> pow(const Real<Option>& s, const Real<Option>& n) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::powf(s.toMachine(), n.toMachine()));
        else {
            static_assert(Option == Float64);
            return Real<Option>(::pow(s.toMachine(), n.toMachine()));
        }
    }
    /*!
     * Ignoring error. If s is a float number, use floor() first.
     * \s must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<ScalarOption Option>
    Real<Option> factorial(const Real<Option>& s) noexcept {
        typedef decltype(s.toMachine()) FloatType;
        const auto trivial = s.toMachine();
        FloatType temp = 1;
        FloatType result = 0;
        while(temp < trivial) {
            temp += 1;
            result *= temp;
        }
        return Real<Option>(result);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> cos(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::cosf(s.toMachine()));
        else
            return Real<Option>(::cos(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> sin(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::sinf(s.toMachine()));
        else
            return Real<Option>(::sin(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline void sincos(Real<Option> s, Real<Option>& sin_result, Real<Option>& cos_result) noexcept {
        using MachineType = Real<Option>::MachineType;
        MachineType sin_temp, cos_temp;
        if constexpr (Option == Double)
            ::sincos(s.toMachine(), (double*)&sin_temp, (double*)&cos_temp);
        else
            ::sincosf(s.toMachine(), (float*)&sin_temp, (float*)&cos_temp);
        sin_result = sin_temp;
        cos_result = cos_temp;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> tan(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::tanf(s.toMachine()));
        else
            return Real<Option>(::tan(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> sec(const Real<Option>& s) noexcept {
        return reciprocal(cos(s));
    }

    template<ScalarOption Option>
    Real<Option> csc(const Real<Option>& s) noexcept {
        return reciprocal(sin(s));
    }

    template<ScalarOption Option>
    Real<Option> cot(const Real<Option>& s) noexcept {
        return cos(s) / sin(s);
    }

    template<ScalarOption Option>
    Real<Option> arccos(const Real<Option>& s) noexcept {
        return Real<Option>(std::acos(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsin(const Real<Option>& s) noexcept {
        return Real<Option>(std::asin(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arctan(const Real<Option>& s) noexcept {
        return Real<Option>(std::atan(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsec(const Real<Option>& s) noexcept {
        return arccos(reciprocal(s));
    }

    template<ScalarOption Option>
    Real<Option> arccsc(const Real<Option>& s) noexcept {
        return arcsin(reciprocal(s));
    }

    template<ScalarOption Option>
    Real<Option> arccot(const Real<Option>& s) noexcept {
        return arctan(reciprocal(s));
    }

    template<ScalarOption Option>
    Real<Option> cosh(const Real<Option>& s) noexcept {
        return Real<Option>(std::cosh(s.toMachine()));
    }


    template<ScalarOption Option>
    Real<Option> sinh(const Real<Option>& s) noexcept {
        return Real<Option>(std::sinh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> tanh(const Real<Option>& s) noexcept {
        return Real<Option>(std::tanh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> sech(const Real<Option>& s) noexcept {
        return Real<Option>(1 / std::cosh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> csch(const Real<Option>& s) noexcept {
        return Real<Option>(1 / std::sinh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> coth(const Real<Option>& s) noexcept {
        return Real<Option>(1 / std::tanh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arccosh(const Real<Option>& s) noexcept {
        return Real<Option>(std::acosh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsinh(const Real<Option>& s) noexcept {
        return Real<Option>(std::asinh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arctanh(const Real<Option>& s) noexcept {
        return Real<Option>(std::atanh(s.toMachine()));
    }

    template<ScalarOption Option>
    Real<Option> arcsech(const Real<Option>& s) noexcept {
        return arccosh(reciprocal(s));
    }

    template<ScalarOption Option>
    Real<Option> arccsch(const Real<Option>& s) noexcept {
        return arcsinh(reciprocal(s));
    }

    template<ScalarOption Option>
    Real<Option> arccoth(const Real<Option>& s) noexcept {
        auto trivial = s.toMachine();
        return Real<Option>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }

    template<ScalarOption Option>
    Real<Option> lncosh(const Real<Option>& s) noexcept {
        using T = Real<Option>;
        const auto x = abs(s);
        return x + ln1p(exp(-x * T(2))) - T(M_LN2);
    }

    template<ScalarOption Option>
    Real<Option> floor(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::floorf(s.toMachine()));
        else
            return Real<Option>(::floor(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Real<Option> ceil(const Real<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Real<Option>(::ceilf(s.toMachine()));
        else
            return Real<Option>(::ceil(s.toMachine()));
    }
}
