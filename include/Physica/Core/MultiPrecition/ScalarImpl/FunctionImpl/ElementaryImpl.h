/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ELEMENTARYIMPL_H
#define PHYSICA_ELEMENTARYIMPL_H

#include "Physica/Core/MultiPrecition/ScalarImpl/ProbabilityFunction.h"
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> square(const Scalar<type, errorTrack>& s) {
        return s * s;
    }
    /*!
     * Calculate @param s * @param s, while this function is faster than simply multiply.
     *
     * Reference: GMP Doc BaseCase Multiplication
     */
    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> square(const Scalar<MultiPrecision, errorTrack>& s) {
        if(s == BasicConst::getInstance().get_1())
            return Scalar(s);
        else {
            const auto s_size = s.getSize();
            //Estimate the length of result first. we will calculate it accurately later.
            auto length = 2 * s_size;
            Scalar<MultiPrecision, errorTrack> result(length, s.power * 2 + 1);

            memset(result.byte, 0, length * sizeof(ScalarUnit));
            for(int i = 0; i < s_size - 1; ++i)
                result[i + s_size] = mulAddArrByWord(result.byte + i + i + 1, s.byte + i + 1, s_size - i - 1, s.byte[i]);
            //Optimize: Shift count is known, possible to optimize the performance.
            //Fix: accuracy is ignored.
            byteLeftShiftEq(result.byte, length, 1);

            ScalarUnit high, low, copy, temp;
            bool carry = false;
            for(unsigned int i = 0; i < s_size; ++i) {
                mulWordByWord(high, low, s.byte[i], s.byte[i]);
                unsigned int double_i = i << 1U;
                /* Handle 2 * i */ {
                    copy = result[double_i];
                    temp = copy + low + carry;
                    carry = copy > temp;
                    result[double_i] = temp;
                }
                /* Handle 2 * i + 1 */ {
                    copy = result[double_i + 1];
                    temp = copy + high + carry;
                    carry = copy > temp;
                    result[double_i + 1] = temp;
                }
            }

            if(high + carry == 0) {
                --result.power;
                --length;
                result.length = length;
                result.byte =
                        reinterpret_cast<ScalarUnit*>(realloc(result.byte, length * sizeof(ScalarUnit)));
            }
            return result;
        }
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> reciprocal(const Scalar<type, errorTrack>& s) {
        return Scalar<type, false>(1) / s;
    }

    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> reciprocal(const Scalar<MultiPrecision, errorTrack>& s) {
        return BasicConst::getInstance().get_1() / s;
    }

    template<ScalarType type>
    Scalar<type, false> sqrt(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::sqrt(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> sqrt(const Scalar<type, true>& s) {
        const auto trivial = s.getTrivial();
        auto error = trivial - s.getA();
        if(error < 0)
            error = trivial + s.getA();
        const auto result = std::sqrt(trivial);
        error = std::sqrt(error) - result;
        error = error < 0 ? -error : error;
        return Scalar<type, true>(result, error);
    }

    template<>
    Scalar<MultiPrecision, false> sqrt(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sqrt(const Scalar<MultiPrecision, true>& s);
    /*!
     * Ignoring error. If s is a float number, use floor() first.
     * \s must be positive, or 1 will be returned.
     *
     * Fix: Easily overflow.
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> factorial(const Scalar<type, errorTrack>& s) {
        const auto trivial = s.getTrivial();
        decltype(trivial) temp = 1;
        decltype(trivial) result;
        while(temp < trivial) {
            temp += 1;
            result *= temp;
        }
        return Scalar<type, errorTrack>(result);
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> factorial(const Scalar<MultiPrecision, errorTrack>& s) {
        //Optimize: Unnecessary copy during floor() if s is a integer itself.
        const Scalar<MultiPrecision, false> integer = floor<false>(s);

        Scalar<MultiPrecision, errorTrack> result(SignedScalarUnit(1));
        Scalar<MultiPrecision, errorTrack> temp(SignedScalarUnit(1));
        while(temp < integer)
            result *= ++temp;
        return result;
    }

    template<ScalarType type>
    Scalar<type, false> ln(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::log(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> ln(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto result = std::log(trivial);
        auto temp = trivial - a;
        auto error = temp < 0 ? result : result - std::log(temp);
        return Scalar<type, true>(result, error);
    }

    template<>
    Scalar<MultiPrecision, false> ln(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> ln(const Scalar<MultiPrecision, true>& s);
    //!Return log_a n
    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    Scalar<type, errorTrack1 | errorTrack2> log(
            const Scalar<type, errorTrack1>& s, const Scalar<type, errorTrack2>& a) {
        return ln(s) / ln(a);
    }

    template<ScalarType type>
    Scalar<type, false> exp(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::exp(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> exp(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto result = std::exp(trivial);
        return Scalar<type, true>(result, std::exp(trivial + s.getA()) - result);
    }

    template<>
    Scalar<MultiPrecision, false> exp(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> exp(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> cos(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::cos(s.getTrivial()));
    }
    //!Accuracy is not accurate.
    template<ScalarType type>
    Scalar<type, true> cos(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto result = std::cos(trivial);
        return Scalar<type, true>(result, std::cos(trivial + s.getA()) - result);
    }

    template<>
    Scalar<MultiPrecision, false> cos(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> cos(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> sin(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::sin(s.getTrivial()));
    }
    //!Accuracy is not accurate.
    template<ScalarType type>
    Scalar<type, true> sin(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto result = std::sin(trivial);
        return Scalar<type, true>(result, std::sin(trivial + s.getA()) - result);
    }

    template<>
    Scalar<MultiPrecision, false> sin(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sin(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> tan(const Scalar<type, errorTrack>& s) {
        return sin(s) / cos(s);
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sec(const Scalar<type, errorTrack>& s) {
        return reciprocal(cos(s));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> csc(const Scalar<type, errorTrack>& s) {
        return reciprocal(sin(s));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> cot(const Scalar<type, errorTrack>& s) {
        return cos(s) / sin(s);
    }

    template<ScalarType type>
    Scalar<type, false> arccos(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::acos(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> arccos(const Scalar<type, true>& s) {
        return Scalar<type, false>(std::acos(s.getTrivial()));
    }
    //!FixIt: arccos() does not consider accuracy.
    template<>
    Scalar<MultiPrecision, false> arccos(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arccos(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> arcsin(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::asin(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> arcsin(const Scalar<type, true>& s) {
        return Scalar<type, false>(std::asin(s.getTrivial()));
    }
    //!FixIt: arcsin() does not consider accuracy.
    template<>
    Scalar<MultiPrecision, false> arcsin(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arcsin(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> arctan(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::atan(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> arctan(const Scalar<type, true>& s) {
        return Scalar<type, false>(std::atan(s.getTrivial()));
    }
    //!FixIt: arctan() does not consider accuracy.
    template<>
    Scalar<MultiPrecision, true> arctan(const Scalar<MultiPrecision, true>& s);

    template<>
    Scalar<MultiPrecision, false> arctan(const Scalar<MultiPrecision, false>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arcsec(const Scalar<type, errorTrack>& s) {
        return arccos(reciprocal(s));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccsc(const Scalar<type, errorTrack>& s) {
        return arcsin(reciprocal(s));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccot(const Scalar<type, errorTrack>& s) {
        return arctan(reciprocal(s));
    }

    template<ScalarType type>
    Scalar<type, false> cosh(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::cosh(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> cosh(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial + a : trivial - a;
        error = std::cosh(error);
        return Scalar<type, true>(std::cosh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> cosh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> cosh(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> sinh(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::sinh(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> sinh(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial + a : trivial - a;
        error = std::sinh(error);
        return Scalar<type, true>(std::sinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> sinh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sinh(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> tanh(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::tanh(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> tanh(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = std::sinh(error);
        return Scalar<type, true>(std::sinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> tanh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> tanh(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> sech(const Scalar<type, false>& s) {
        return Scalar<type, false>(1 / std::cosh(s.getTrivial()));
    }

    template<ScalarType type>
    Scalar<type, true> sech(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = 1 / std::cosh(error);
        return Scalar<type, true>(1 / std::cosh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> sech(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> sech(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> csch(const Scalar<type, false>& s) {
        return Scalar<type, false>(1 / std::sinh(s.getTrivial()));
    }
    //!FixIt: Accuracy is not accurate if trivial <= a.
    template<ScalarType type>
    Scalar<type, true> csch(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = 1 / std::sinh(error);
        return Scalar<type, true>(1 / std::sinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> csch(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> csch(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> coth(const Scalar<type, false>& s) {
        return Scalar<type, false>(1 / std::tanh(s.getTrivial()));
    }
    //!FixIt: Accuracy is not accurate if trivial <= a.
    template<ScalarType type>
    Scalar<type, true> coth(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial < 0 ? trivial + a : trivial - a;
        error = 1 / std::tanh(error);
        return Scalar<type, true>(1 / std::tanh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> coth(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> coth(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> arccosh(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::acosh(s.getTrivial()));
    }
    //!Accuracy is not accurate when trivial <= a.
    template<ScalarType type>
    Scalar<type, true> arccosh(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > a ? trivial - a : trivial + a;
        error = acosh(error);
        return Scalar<type, true>(std::acosh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> arccosh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arccosh(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> arcsinh(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::asinh(s.getTrivial()));
    }
    //!Accuracy is not accurate when trivial <= a.
    template<ScalarType type>
    Scalar<type, true> arcsinh(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial - a : trivial + a;
        error = asinh(error);
        return Scalar<type, true>(std::asinh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> arcsinh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arcsinh(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type>
    Scalar<type, false> arctanh(const Scalar<type, false>& s) {
        return Scalar<type, false>(std::atanh(s.getTrivial()));
    }
    //!Accuracy is not accurate when trivial + a > 1 or trivial - a < -1.
    template<ScalarType type>
    Scalar<type, true> arctanh(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial + a : trivial - a;
        error = atanh(error);
        return Scalar<type, true>(std::atanh(trivial), error);
    }

    template<>
    Scalar<MultiPrecision, false> arctanh(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arctanh(const Scalar<MultiPrecision, true>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arcsech(const Scalar<type, errorTrack>& s) {
        return arccosh(reciprocal(s));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arccsch(const Scalar<type, errorTrack>& s) {
        return arcsinh(reciprocal(s));
    }

    template<ScalarType type>
    Scalar<type, false> arccoth(const Scalar<type, false>& s) {
        auto trivial = s.getTrivial();
        return Scalar<type, false>(std::log((trivial + 1) / (trivial - 1)) / 2);
    }
    //!Accuracy is not accurate when |trivial| < a.
    template<ScalarType type>
    Scalar<type, true> arccoth(const Scalar<type, true>& s) {
        auto trivial = s.getTrivial();
        auto a = s.getA();
        auto error = trivial > 0 ? trivial - a : trivial + a;
        error = std::log((error + 1) / (error - 1)) / 2;
        return Scalar<type, false>(std::log((trivial + 1) / (trivial - 1)) / 2, error);
    }

    template<>
    Scalar<MultiPrecision, false> arccoth(const Scalar<MultiPrecision, false>& s);

    template<>
    Scalar<MultiPrecision, true> arccoth(const Scalar<MultiPrecision, true>& s);
}

#endif