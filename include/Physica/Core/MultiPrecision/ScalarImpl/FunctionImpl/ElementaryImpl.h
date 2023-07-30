/*
 * Copyright 2020-2023 WeiBo He.
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

#include "Physica/Core/MultiPrecision/ScalarImpl/ProbabilityFunction.h"
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<ScalarOption option>
    __host__ __device__ inline Scalar<option> abs(const Scalar<option>& s) {
        using T = Scalar<option>;
        T temp(s);
        return T(std::move(temp.toAbs()));
    }

    template<ScalarOption option>
    __host__ __device__ Scalar<option> square(const Scalar<option>& s) {
        return s * s;
    }
    /**
     * Calculate @param s * @param s, while this function is faster than simply multiply.
     *
     * Reference: GMP Doc BaseCase Multiplication
     */
    template<>
    Scalar<MultiPrecision> square(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    __host__ __device__ inline Scalar<option> reciprocal(const Scalar<option>& s) {
        return Scalar<option>(1) / s;
    }

    template<>
    inline Scalar<MultiPrecision> reciprocal(const Scalar<MultiPrecision>& s) {
        return BasicConst::getInstance()._1 / s;
    }

    template<ScalarOption option>
    __host__ __device__ Scalar<option> sqrt(const Scalar<option>& s) {
        return Scalar<option>(std::sqrt(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option> cbrt(const Scalar<option>& s) {
        return Scalar<option>(std::cbrt(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> sqrt(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> ln(const Scalar<option>& s) {
        return Scalar<option>(std::log(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> ln(const Scalar<MultiPrecision>& s);
    /**
     * \return log_a n
     */
    template<ScalarOption option>
    Scalar<option> log(const Scalar<option>& s, const Scalar<option>& a) {
        return ln(s) / ln(a);
    }

    template<ScalarOption option>
    Scalar<option> exp(const Scalar<option>& s) {
        return Scalar<option>(std::exp(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> exp(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> pow(const Scalar<option>& s, const Scalar<option>& n) {
        if constexpr (option == MultiPrecision)
            return exp(n * ln(s));
        else
            return Scalar<option>(std::pow(s.getTrivial(), n.getTrivial()));
    }
    /*!
     * Ignoring error. If s is a float number, use floor() first.
     * \s must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<ScalarOption option>
    Scalar<option> factorial(const Scalar<option>& s) {
        if constexpr (option == MultiPrecision) {
            //Optimize: Unnecessary copy during floor() if s is a integer itself.
            const Scalar<MultiPrecision> integer = floor(s);

            Scalar<MultiPrecision> result(SignedMPUnit(1));
            Scalar<MultiPrecision> temp(SignedMPUnit(1));
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
            return Scalar<option>(result);
        }
    }

    template<>
    Scalar<MultiPrecision> factorial(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> cos(const Scalar<option>& s) {
        return Scalar<option>(std::cos(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> cos(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> sin(const Scalar<option>& s) {
        return Scalar<option>(std::sin(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> sin(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    void sincos(Scalar<option> s, Scalar<option>& sin_result, Scalar<option>& cos_result) {
        using TrivialType = typename Scalar<option>::TrivialType;
        TrivialType sin_temp, cos_temp;
        if constexpr (option == Double)
            ::sincos(s.getTrivial(), (double*)&sin_temp, (double*)&cos_temp);
        else
            ::sincosf(s.getTrivial(), (float*)&sin_temp, (float*)&cos_temp);
        sin_result = sin_temp;
        cos_result = cos_temp;
    }

    template<ScalarOption option>
    Scalar<option> tan(const Scalar<option>& s) {
        return sin(s) / cos(s);
    }

    template<ScalarOption option>
    Scalar<option> sec(const Scalar<option>& s) {
        return reciprocal(cos(s));
    }

    template<ScalarOption option>
    Scalar<option> csc(const Scalar<option>& s) {
        return reciprocal(sin(s));
    }

    template<ScalarOption option>
    Scalar<option> cot(const Scalar<option>& s) {
        return cos(s) / sin(s);
    }

    template<ScalarOption option>
    Scalar<option> arccos(const Scalar<option>& s) {
        return Scalar<option>(std::acos(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> arccos(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> arcsin(const Scalar<option>& s) {
        return Scalar<option>(std::asin(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> arcsin(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> arctan(const Scalar<option>& s) {
        return Scalar<option>(std::atan(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> arctan(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> arcsec(const Scalar<option>& s) {
        return arccos(reciprocal(s));
    }

    template<ScalarOption option>
    Scalar<option> arccsc(const Scalar<option>& s) {
        return arcsin(reciprocal(s));
    }

    template<ScalarOption option>
    Scalar<option> arccot(const Scalar<option>& s) {
        return arctan(reciprocal(s));
    }

    template<ScalarOption option>
    Scalar<option> cosh(const Scalar<option>& s) {
        return Scalar<option>(std::cosh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> cosh(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> sinh(const Scalar<option>& s) {
        return Scalar<option>(std::sinh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> sinh(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> tanh(const Scalar<option>& s) {
        return Scalar<option>(std::tanh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> tanh(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> sech(const Scalar<option>& s) {
        return Scalar<option>(1 / std::cosh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> sech(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> csch(const Scalar<option>& s) {
        return Scalar<option>(1 / std::sinh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> csch(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> coth(const Scalar<option>& s) {
        return Scalar<option>(1 / std::tanh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> coth(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> arccosh(const Scalar<option>& s) {
        return Scalar<option>(std::acosh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> arccosh(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> arcsinh(const Scalar<option>& s) {
        return Scalar<option>(std::asinh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> arcsinh(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> arctanh(const Scalar<option>& s) {
        return Scalar<option>(std::atanh(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> arctanh(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    Scalar<option> arcsech(const Scalar<option>& s) {
        return arccosh(reciprocal(s));
    }

    template<ScalarOption option>
    Scalar<option> arccsch(const Scalar<option>& s) {
        return arcsinh(reciprocal(s));
    }

    template<ScalarOption option>
    Scalar<option> arccoth(const Scalar<option>& s) {
        auto trivial = s.getTrivial();
        return Scalar<option>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }

    template<>
    Scalar<MultiPrecision> arccoth(const Scalar<MultiPrecision>& s);
}
