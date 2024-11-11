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
    __host__ __device__ inline Scalar<Option> abs(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::fabsf(s.toMachine()));
        else
            return Scalar<Option>(::fabs(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> relu(const Scalar<Option>& s) noexcept {
        return s.isPositive() ? s : Scalar<Option>(0);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> square(const Scalar<Option>& s) noexcept {
        return s * s;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> reciprocal(const Scalar<Option>& s) noexcept {
        assert(!s.isZero() && "[Error]: Divide by zero");
        return Scalar<Option>(1) / s;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> sqrt(const Scalar<Option>& s) noexcept {
        assert(!s.isNegative());
        if constexpr (Option == Float)
            return Scalar<Option>(::sqrtf(s.toMachine()));
        else
            return Scalar<Option>(::sqrt(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> cbrt(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::cbrtf(s.toMachine()));
        else
            return Scalar<Option>(::cbrt(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ln(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::logf(s.toMachine()));
        else
            return Scalar<Option>(::log(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ln1p(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::log1pf(s.toMachine()));
        else
            return Scalar<Option>(::log1p(s.toMachine()));
    }
    /**
     * \return log_a n
     */
    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> log(const Scalar<Option>& s, const Scalar<Option>& a) noexcept {
        return ln(s) / ln(a);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> exp(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::expf(s.toMachine()));
        else
            return Scalar<Option>(::exp(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> pow(const Scalar<Option>& s, const Scalar<Option>& n) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::powf(s.toMachine(), n.toMachine()));
        else {
            static_assert(Option == Float64);
            return Scalar<Option>(::pow(s.toMachine(), n.toMachine()));
        }
    }
    /*!
     * Ignoring error. If s is a float number, use floor() first.
     * \s must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<ScalarOption Option>
    Scalar<Option> factorial(const Scalar<Option>& s) noexcept {
        typedef decltype(s.toMachine()) FloatType;
        const auto trivial = s.toMachine();
        FloatType temp = 1;
        FloatType result = 0;
        while(temp < trivial) {
            temp += 1;
            result *= temp;
        }
        return Scalar<Option>(result);
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> cos(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::cosf(s.toMachine()));
        else
            return Scalar<Option>(::cos(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> sin(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::sinf(s.toMachine()));
        else
            return Scalar<Option>(::sin(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline void sincos(Scalar<Option> s, Scalar<Option>& sin_result, Scalar<Option>& cos_result) noexcept {
        using MachineType = typename Scalar<Option>::MachineType;
        MachineType sin_temp, cos_temp;
        if constexpr (Option == Double)
            ::sincos(s.toMachine(), (double*)&sin_temp, (double*)&cos_temp);
        else
            ::sincosf(s.toMachine(), (float*)&sin_temp, (float*)&cos_temp);
        sin_result = sin_temp;
        cos_result = cos_temp;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> tan(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::tanf(s.toMachine()));
        else
            return Scalar<Option>(::tan(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> sec(const Scalar<Option>& s) noexcept {
        return reciprocal(cos(s));
    }

    template<ScalarOption Option>
    Scalar<Option> csc(const Scalar<Option>& s) noexcept {
        return reciprocal(sin(s));
    }

    template<ScalarOption Option>
    Scalar<Option> cot(const Scalar<Option>& s) noexcept {
        return cos(s) / sin(s);
    }

    template<ScalarOption Option>
    Scalar<Option> arccos(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::acos(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> arcsin(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::asin(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> arctan(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::atan(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> arcsec(const Scalar<Option>& s) noexcept {
        return arccos(reciprocal(s));
    }

    template<ScalarOption Option>
    Scalar<Option> arccsc(const Scalar<Option>& s) noexcept {
        return arcsin(reciprocal(s));
    }

    template<ScalarOption Option>
    Scalar<Option> arccot(const Scalar<Option>& s) noexcept {
        return arctan(reciprocal(s));
    }

    template<ScalarOption Option>
    Scalar<Option> cosh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::cosh(s.toMachine()));
    }


    template<ScalarOption Option>
    Scalar<Option> sinh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::sinh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> tanh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::tanh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> sech(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(1 / std::cosh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> csch(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(1 / std::sinh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> coth(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(1 / std::tanh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> arccosh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::acosh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> arcsinh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::asinh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> arctanh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::atanh(s.toMachine()));
    }

    template<ScalarOption Option>
    Scalar<Option> arcsech(const Scalar<Option>& s) noexcept {
        return arccosh(reciprocal(s));
    }

    template<ScalarOption Option>
    Scalar<Option> arccsch(const Scalar<Option>& s) noexcept {
        return arcsinh(reciprocal(s));
    }

    template<ScalarOption Option>
    Scalar<Option> arccoth(const Scalar<Option>& s) noexcept {
        auto trivial = s.toMachine();
        return Scalar<Option>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }

    template<ScalarOption Option>
    Scalar<Option> lncosh(const Scalar<Option>& s) noexcept {
        using T = Scalar<Option>;
        const auto x = abs(s);
        return x + ln1p(exp(-x * T(2))) - T(M_LN2);
    }

    template<ScalarOption Option>
    Scalar<Option> floor(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::floorf(s.toMachine()));
        else
            return Scalar<Option>(::floor(s.toMachine()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ceil(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::ceilf(s.toMachine()));
        else
            return Scalar<Option>(::ceil(s.toMachine()));
    }
}
