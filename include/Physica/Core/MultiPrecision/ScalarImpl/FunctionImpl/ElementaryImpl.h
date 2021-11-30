/*
 * Copyright 2020-2021 WeiBo He.
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
#ifndef PHYSICA_ELEMENTARYIMPL_H
#define PHYSICA_ELEMENTARYIMPL_H

#include "Physica/Core/MultiPrecision/ScalarImpl/ProbabilityFunction.h"
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> abs(const Scalar<option, errorTrack>& s) {
        using T = Scalar<option, errorTrack>;
        T temp(s);
        return T(std::move(temp.toAbs()));
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> square(const Scalar<option, errorTrack>& s) {
        return s * s;
    }
    /*!
     * Calculate @param s * @param s, while this function is faster than simply multiply.
     *
     * Reference: GMP Doc BaseCase Multiplication
     */
    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> square(const Scalar<MultiPrecision, errorTrack>& s) {
        if(s == BasicConst::getInstance()._1)
            return Scalar(s);
        else {
            const auto s_size = s.getSize();
            //Estimate the length of result first. we will calculate it accurately later.
            auto length = 2 * s_size;
            Scalar<MultiPrecision, errorTrack> result(length, s.getPower() * 2 + 1);

            memset(result.byte, 0, length * sizeof(MPUnit));
            for(int i = 0; i < s_size - 1; ++i)
                result.setByte(i + s_size
                               , mulAddArrByWord(result.byte + i + i + 1, s.byte + i + 1, s_size - i - 1, s.byte[i]));
            //Optimize: Shift count is known, possible to optimize the performance.
            //Fix: accuracy is ignored.
            byteLeftShiftEq(result.byte, length, 1);

            MPUnit high{}, low, copy, temp;
            bool carry = false;
            for(int i = 0; i < s_size; ++i) {
                mulWordByWord(high, low, s.byte[i], s.byte[i]);
                unsigned int double_i = static_cast<unsigned int>(i) << 1U;
                /* Handle 2 * i */ {
                    copy = result[double_i];
                    temp = copy + low + carry;
                    carry = copy > temp;
                    result.setByte(double_i, temp);
                }
                /* Handle 2 * i + 1 */ {
                    copy = result[double_i + 1];
                    temp = copy + high + carry;
                    carry = copy > temp;
                    result.setByte(double_i + 1, temp);
                }
            }

            if(high + carry == 0) {
                --result.power;
                --length;
                result.length = length;
                result.byte =
                        reinterpret_cast<MPUnit*>(realloc(result.byte, length * sizeof(MPUnit)));
            }
            return result;
        }
    }

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> reciprocal(const Scalar<option, errorTrack>& s) {
        return Scalar<option, false>(1) / s;
    }

    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> reciprocal(const Scalar<MultiPrecision, errorTrack>& s) {
        return BasicConst::getInstance()._1 / s;
    }

    template<ScalarOption option>
    Scalar<option, false> sqrt(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::sqrt(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> sqrt(const Scalar<option, true>& s) {
        const auto trivial = s.getTrivial();
        auto error = trivial - s.getA();
        if(error < 0)
            error = trivial + s.getA();
        const auto result = std::sqrt(trivial);
        error = std::sqrt(error) - result;
        error = error < 0 ? -error : error;
        return Scalar<option, true>(result, error);
    }

    template<ScalarOption option>
    Scalar<option, false> cbrt(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::cbrt(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision, false> sqrt(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sqrt(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> pow(const Scalar<option, errorTrack>& s, const Scalar<option, errorTrack>& n) {
        return Scalar<option, errorTrack>(std::pow(s.getTrivial(), n.getTrivial()));
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> pow(const Scalar<MultiPrecision, errorTrack>& s, const Scalar<MultiPrecision, errorTrack>& n) {
        return exp(n * ln(s));
    }
    /*!
     * Ignoring error. If s is a float number, use floor() first.
     * \s must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> factorial(const Scalar<option, errorTrack>& s) {
        typedef decltype(s.getTrivial()) FloatType;
        const auto trivial = s.getTrivial();
        FloatType temp = 1;
        FloatType result = 0;
        while(temp < trivial) {
            temp += 1;
            result *= temp;
        }
        return Scalar<option, errorTrack>(result);
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> factorial(const Scalar<MultiPrecision, errorTrack>& s) {
        //Optimize: Unnecessary copy during floor() if s is a integer itself.
        const Scalar<MultiPrecision, false> integer = floor<false>(s);

        Scalar<MultiPrecision, errorTrack> result(SignedMPUnit(1));
        Scalar<MultiPrecision, errorTrack> temp(SignedMPUnit(1));
        while(temp < integer)
            result *= ++temp;
        return result;
    }

    template<ScalarOption option>
    Scalar<option, false> ln(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::log(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> ln(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto result = std::log(trivial);
        auto temp = trivial - a;
        auto error = temp < 0 ? result : result - std::log(temp);
        return Scalar<option, true>(result, error);
    }

    template<>
    Scalar<MultiPrecision, false> ln(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> ln(const Scalar<MultiPrecision, true>& s);
    //!Return log_a n
    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    Scalar<option, errorTrack1 || errorTrack2> log(
            const Scalar<option, errorTrack1>& s, const Scalar<option, errorTrack2>& a) {
        return ln(s) / ln(a);
    }

    template<ScalarOption option>
    Scalar<option, false> exp(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::exp(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> exp(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto result = std::exp(trivial);
        return Scalar<option, true>(result, std::exp(trivial + s.getA()) - result);
    }

    template<>
    Scalar<MultiPrecision, false> exp(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> exp(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> cos(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::cos(s.getTrivial()));
    }
    //!Accuracy is not accurate.
    template<ScalarOption option>
    Scalar<option, true> cos(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto result = std::cos(trivial);
        return Scalar<option, true>(result, std::cos(trivial + s.getA()) - result);
    }

    template<>
    Scalar<MultiPrecision, false> cos(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> cos(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> sin(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::sin(s.getTrivial()));
    }
    //!Accuracy is not accurate.
    template<ScalarOption option>
    Scalar<option, true> sin(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto result = std::sin(trivial);
        return Scalar<option, true>(result, std::sin(trivial + s.getA()) - result);
    }

    template<>
    Scalar<MultiPrecision, false> sin(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sin(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> tan(const Scalar<option, errorTrack>& s) {
        return sin(s) / cos(s);
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sec(const Scalar<option, errorTrack>& s) {
        return reciprocal(cos(s));
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> csc(const Scalar<option, errorTrack>& s) {
        return reciprocal(sin(s));
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> cot(const Scalar<option, errorTrack>& s) {
        return cos(s) / sin(s);
    }

    template<ScalarOption option>
    Scalar<option, false> arccos(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::acos(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> arccos(const Scalar<option, true>& s) {
        return Scalar<option, false>(std::acos(s.getTrivial()));
    }
    //!FixIt: arccos() does not consider accuracy.
    template<>
    Scalar<MultiPrecision, false> arccos(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arccos(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> arcsin(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::asin(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> arcsin(const Scalar<option, true>& s) {
        return Scalar<option, false>(std::asin(s.getTrivial()));
    }
    //!FixIt: arcsin() does not consider accuracy.
    template<>
    Scalar<MultiPrecision, false> arcsin(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arcsin(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> arctan(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::atan(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> arctan(const Scalar<option, true>& s) {
        return Scalar<option, false>(std::atan(s.getTrivial()));
    }
    //!FixIt: arctan() does not consider accuracy.
    template<>
    Scalar<MultiPrecision, true> arctan(const Scalar<MultiPrecision, true>& s);

    template<>
    Scalar<MultiPrecision, false> arctan(const Scalar<MultiPrecision, false>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arcsec(const Scalar<option, errorTrack>& s) {
        return arccos(reciprocal(s));
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccsc(const Scalar<option, errorTrack>& s) {
        return arcsin(reciprocal(s));
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccot(const Scalar<option, errorTrack>& s) {
        return arctan(reciprocal(s));
    }

    template<ScalarOption option>
    Scalar<option, false> cosh(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::cosh(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> cosh(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial + a : trivial - a;
        error = std::cosh(error);
        return Scalar<option, true>(std::cosh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> cosh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> cosh(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> sinh(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::sinh(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> sinh(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial + a : trivial - a;
        error = std::sinh(error);
        return Scalar<option, true>(std::sinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> sinh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sinh(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> tanh(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::tanh(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> tanh(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = std::sinh(error);
        return Scalar<option, true>(std::sinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> tanh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> tanh(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> sech(const Scalar<option, false>& s) {
        return Scalar<option, false>(1 / std::cosh(s.getTrivial()));
    }

    template<ScalarOption option>
    Scalar<option, true> sech(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = 1 / std::cosh(error);
        return Scalar<option, true>(1 / std::cosh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> sech(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sech(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> csch(const Scalar<option, false>& s) {
        return Scalar<option, false>(1 / std::sinh(s.getTrivial()));
    }
    //!FixIt: Accuracy is not accurate if trivial <= a.
    template<ScalarOption option>
    Scalar<option, true> csch(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = 1 / std::sinh(error);
        return Scalar<option, true>(1 / std::sinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> csch(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> csch(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> coth(const Scalar<option, false>& s) {
        return Scalar<option, false>(1 / std::tanh(s.getTrivial()));
    }
    //!FixIt: Accuracy is not accurate if trivial <= a.
    template<ScalarOption option>
    Scalar<option, true> coth(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = 1 / std::tanh(error);
        return Scalar<option, true>(1 / std::tanh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> coth(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> coth(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> arccosh(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::acosh(s.getTrivial()));
    }
    //!Accuracy is not accurate when trivial <= a.
    template<ScalarOption option>
    Scalar<option, true> arccosh(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > a ? trivial - a : trivial + a;
        error = acosh(error);
        return Scalar<option, true>(std::acosh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> arccosh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arccosh(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> arcsinh(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::asinh(s.getTrivial()));
    }
    //!Accuracy is not accurate when trivial <= a.
    template<ScalarOption option>
    Scalar<option, true> arcsinh(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial - a : trivial + a;
        error = asinh(error);
        return Scalar<option, true>(std::asinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> arcsinh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arcsinh(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option>
    Scalar<option, false> arctanh(const Scalar<option, false>& s) {
        return Scalar<option, false>(std::atanh(s.getTrivial()));
    }
    //!Accuracy is not accurate when trivial + a > 1 or trivial - a < -1.
    template<ScalarOption option>
    Scalar<option, true> arctanh(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial + a : trivial - a;
        error = atanh(error);
        return Scalar<option, true>(std::atanh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> arctanh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arctanh(const Scalar<MultiPrecision, true>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arcsech(const Scalar<option, errorTrack>& s) {
        return arccosh(reciprocal(s));
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arccsch(const Scalar<option, errorTrack>& s) {
        return arcsinh(reciprocal(s));
    }

    template<ScalarOption option>
    Scalar<option, false> arccoth(const Scalar<option, false>& s) {
        auto trivial = s.getTrivial();
        return Scalar<option, false>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }
    //!Accuracy is not accurate when |trivial| < a.
    template<ScalarOption option>
    Scalar<option, true> arccoth(const Scalar<option, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial - a : trivial + a;
        error = std::log((error + 1) / (error - 1)) / 2;
        return Scalar<option, false>(std::log((trivial + 1) / (trivial - 1)) / 2, error);
    }

    template<>
    Scalar<MultiPrecision, false> arccoth(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arccoth(const Scalar<MultiPrecision, true>& s);
}

#endif