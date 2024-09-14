/*
 * Copyright 2020-2024 Weibo He.
 *
 * This file is part of Physica.

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

#include <Physica/Core/MultiPrecision/ScalarImpl/ProbabilityFunction.h>
#include <Physica/Core/MultiPrecision/ScalarImpl/Operation/Pow.h>

namespace Physica::Core {
    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> abs(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::fabsf(s.getTrivial()));
        else
            return Scalar<Option>(::fabs(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> abs(const Scalar<FloatMP>& s) noexcept;

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
        return Scalar<Option>(1) / s;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> sqrt(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::sqrtf(s.getTrivial()));
        else
            return Scalar<Option>(::sqrt(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> sqrt(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> cbrt(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::cbrtf(s.getTrivial()));
        else
            return Scalar<Option>(::cbrt(s.getTrivial()));
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ln(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::logf(s.getTrivial()));
        else
            return Scalar<Option>(::log(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> ln(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> ln1p(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::log1pf(s.getTrivial()));
        else
            return Scalar<Option>(::log1p(s.getTrivial()));
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
            return Scalar<Option>(::expf(s.getTrivial()));
        else
            return Scalar<Option>(::exp(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> exp(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> pow(const Scalar<Option>& s, const Scalar<Option>& n) noexcept {
        if constexpr (Option == FloatMP)
            return exp(n * ln(s));
        else if constexpr (Option == Float)
            return Scalar<Option>(::powf(s.getTrivial(), n.getTrivial()));
        else
            return Scalar<Option>(::pow(s.getTrivial(), n.getTrivial()));
    }

    template<>
    inline Scalar<FloatMP> pow(const Scalar<FloatMP>& s1, const Scalar<FloatMP>& s2) noexcept {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }
    /*!
     * Ignoring error. If s is a float number, use floor() first.
     * \s must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<ScalarOption Option>
    Scalar<Option> factorial(const Scalar<Option>& s) noexcept {
        if constexpr (Option == FloatMP) {
            //Optimize: Unnecessary copy during floor() if s is a integer itself.
            const Scalar<FloatMP> integer = floor(s);

            Scalar<FloatMP> result(SignedMPUnit(1));
            Scalar<FloatMP> temp(SignedMPUnit(1));
            while(temp < integer)
                result *= ++temp;
            return result;
        }
        else {
            typedef decltype(s.getTrivial()) FloatType;
            const auto trivial = s.getTrivial();
            FloatType temp = 1;
            FloatType result = 0;
            while(temp < trivial) {
                temp += 1;
                result *= temp;
            }
            return Scalar<Option>(result);
        }
    }

    template<>
    Scalar<FloatMP> factorial(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> cos(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::cosf(s.getTrivial()));
        else
            return Scalar<Option>(::cos(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> cos(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> sin(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::sinf(s.getTrivial()));
        else
            return Scalar<Option>(::sin(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> sin(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    __host__ __device__ inline void sincos(Scalar<Option> s, Scalar<Option>& sin_result, Scalar<Option>& cos_result) noexcept {
        using TrivialType = typename Scalar<Option>::TrivialType;
        TrivialType sin_temp, cos_temp;
        if constexpr (Option == Double)
            ::sincos(s.getTrivial(), (double*)&sin_temp, (double*)&cos_temp);
        else
            ::sincosf(s.getTrivial(), (float*)&sin_temp, (float*)&cos_temp);
        sin_result = sin_temp;
        cos_result = cos_temp;
    }

    template<ScalarOption Option>
    __host__ __device__ inline Scalar<Option> tan(const Scalar<Option>& s) noexcept {
        if constexpr (Option == Float)
            return Scalar<Option>(::tanf(s.getTrivial()));
        else
            return Scalar<Option>(::tan(s.getTrivial()));
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
        return Scalar<Option>(std::acos(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> arccos(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arcsin(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::asin(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> arcsin(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arctan(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::atan(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> arctan(const Scalar<FloatMP>& s) noexcept;

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
        return Scalar<Option>(std::cosh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> cosh(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> sinh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::sinh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> sinh(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> tanh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::tanh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> tanh(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> sech(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(1 / std::cosh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> sech(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> csch(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(1 / std::sinh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> csch(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> coth(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(1 / std::tanh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> coth(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arccosh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::acosh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> arccosh(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arcsinh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::asinh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> arcsinh(const Scalar<FloatMP>& s) noexcept;

    template<ScalarOption Option>
    Scalar<Option> arctanh(const Scalar<Option>& s) noexcept {
        return Scalar<Option>(std::atanh(s.getTrivial()));
    }

    template<>
    Scalar<FloatMP> arctanh(const Scalar<FloatMP>& s) noexcept;

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
        auto trivial = s.getTrivial();
        return Scalar<Option>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }

    template<>
    Scalar<FloatMP> arccoth(const Scalar<FloatMP>& s) noexcept;
}
