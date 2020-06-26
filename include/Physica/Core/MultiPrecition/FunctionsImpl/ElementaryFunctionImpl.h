/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ELEMENTARYFUNCTIONIMPL_H
#define PHYSICA_ELEMENTARYFUNCTIONIMPL_H

#include <cmath>
#include "Physica/Core/Math/Algebra/BasicAlgebra/Solve.h"
#include "Physica/Core/MultiPrecition/ProbabilityFunction.h"
/*!
 * This file is part of implementations of \ElementaryFunction.
 * Do not include this header file, include ElementaryFunction.h instead.
 */
namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> square(const Scalar<type, errorTrack>& s) {
        return s * s;
    }
    //!Reference: GMP Doc BaseCase Multiplication
    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> square(const Scalar<MultiPrecision, errorTrack>& s) {
        if(s == BasicConst::getInstance().get_1())
            return Scalar(s);
        else {
            auto s_size = s.getSize();
            //Estimate the ed of result first. we will calculate it accurately later.
            const auto length = 2 * s_size;
            Scalar<MultiPrecision, errorTrack> result(length, s.power * 2 + 1);

            for(int i = 0; i < s_size - 1; ++i)
                result[i + s_size] = mulAddArrByWord(result.byte + i + i + 1
                        , s.byte + i + 1, s_size - i - 1, s.byte[i]);
            //Optimize: Shift count is known, possible to optimize the performance.
            //Fix: accuracy is ignored.
            byteLeftShiftEq(result.byte, length, 1);

            ScalarUnit high, low;
            for(int i = 0; i < s_size; ++i) {
                mulWordByWord(high, low, s.byte[i], s.byte[i]);
                result[2 * i] += low;
                result[2 * i + 1] += high;
            }

            if(high == 0) {
                --result.power;
                result.byte =
                        reinterpret_cast<ScalarUnit*>(realloc(result.byte, (length - 1) * sizeof(ScalarUnit)));
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
    Scalar<MultiPrecision, false> sqrt(const Scalar<MultiPrecision, false>& s) {
        if(s.isNegative())
            qFatal("Can not resolve the square root of a minus value.");
        if(s.isZero())
            return Scalar(BasicConst::getInstance().get_0());
        Scalar copy_s(s);
        //Let s < 1 so as to control error.
        int add_power = 0;
        if(copy_s.getPower() > 0) {
            if(copy_s.getPower() % 2 == 0) {
                add_power = copy_s.getPower() / 2 + 1;
                copy_s.power = -2;
            }
            else {
                add_power = (copy_s.getPower() + 1) / 2;
                copy_s.power = -1;
            }
        }

        Scalar result = Scalar<MultiPrecision, false>(static_cast<SignedScalarUnit>(1));
        //3.33 is the big approximate value of ln(10)/ln(2)
        for(int i = 0; i < LONG_WIDTH * GlobalPrecision; ++i)
            result = (result + Scalar<MultiPrecision, false>::div<false>(copy_s, result)) >> 1U;
        result.power += add_power;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> sqrt(const Scalar<MultiPrecision, true>& s) {
        if(s.isNegative())
            qFatal("Can not resolve the square root of a minus value.");
        if(s.isZero())
            return Scalar<MultiPrecision, true>(BasicConst::getInstance().get_0());
        Scalar copy_s(s);
        //Let s < 1 so as to control error.
        int add_power = 0;
        if(copy_s.getPower() > 0) {
            if(copy_s.getPower() % 2 == 0) {
                add_power = copy_s.getPower() / 2 + 1;
                copy_s.power = -2;
            }
            else {
                add_power = (copy_s.getPower() + 1) / 2;
                copy_s.power = -1;
            }
        }

        Scalar result = Scalar<MultiPrecision, true>(static_cast<SignedScalarUnit>(1));
        //3.33 is the big approximate value of ln(10)/ln(2)
        for(int i = 0; i < LONG_WIDTH * GlobalPrecision; ++i)
            result = (result + Scalar<MultiPrecision, false>::div<false>(copy_s, result)) >> 1U;
        result.power += add_power;

        result.toUnitA();
        if(s.getA() != 0) {
            Scalar<MultiPrecision, false> s_error = s.getMinimum();
            if(s_error.isNegative())
                s_error = s.getMaximum();
            Scalar<MultiPrecision, false> error = sqrt(s_error);
            error -= result;
            error.toAbs();

            result.applyError(error);
        }
        return result;
    }
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
        const auto& _1 = BasicConst::getInstance().get_1();

        Scalar<MultiPrecision, errorTrack> result(SignedScalarUnit(1));
        Scalar<MultiPrecision, errorTrack> temp(SignedScalarUnit(1));
        while(temp < integer) {
            temp += _1;
            result *= temp;
        }
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
    Scalar<MultiPrecision, false> ln(const Scalar<MultiPrecision, false>& s) {
        if(!s.isPositive())
            qFatal("Can not resolve the logarithm of zero or a negative value.");
        Scalar<MultiPrecision, false> result(static_cast<SignedScalarUnit>(0));
        if(s == BasicConst::getInstance().get_1())
            return result;
        const auto& _1 = BasicConst::getInstance().get_1();
        auto temp_1 = Scalar<MultiPrecision, false>::sub<false>(s, _1)
                      / Scalar<MultiPrecision, false>::add<false>(s, _1);
        Scalar<MultiPrecision, false> copy_temp_1(temp_1);
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(1));

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / rank;
            result += temp;

            temp_1 *= copy_temp_1;
            rank += BasicConst::getInstance().get_1();
            Scalar criteria = temp_1 / rank;
            //Break if result meets the precision goal.
            if(result.getPower() - criteria.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= copy_temp_1;
            rank += _1;
        }
        result *= BasicConst::getInstance().get_2();
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> ln(const Scalar<MultiPrecision, true>& s) {
        if(!s.isPositive())
            qFatal("Can not resolve the logarithm of zero or a negative value.");
        Scalar<MultiPrecision, true> result(static_cast<SignedScalarUnit>(0));
        if(s == BasicConst::getInstance().get_1())
            return result;
        const auto& _1 = BasicConst::getInstance().get_1();
        auto temp_1 = Scalar<MultiPrecision, false>::sub<true>(s, _1)
                      / Scalar<MultiPrecision, false>::add<false>(s, _1);
        Scalar<MultiPrecision, true> copy_temp_1(temp_1);
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(1));

        temp_1.toUnitA();
        copy_temp_1.toUnitA();
        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / rank;
            temp.clearA();
            result += temp;

            temp_1 *= copy_temp_1;
            rank += BasicConst::getInstance().get_1();
            Scalar criteria = temp_1 / rank;
            //Break if result meets the precision goal.
            if(result.getPower() - criteria.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= copy_temp_1;
            rank += _1;
        }
        result *= BasicConst::getInstance().get_2();

        if(s.getA() != 0) {
            Scalar<MultiPrecision, false> s_error = s.getMinimum();
            if(s_error.isNegative())
                s_error = s.getMaximum();
            Scalar<MultiPrecision, false> error = ln(s_error);
            error -= result;
            error.toAbs();

            result.applyError(error);
        }
        return result;
    }
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
    Scalar<MultiPrecision, false> exp(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedScalarUnit>(1));
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(1));
        Scalar<MultiPrecision, false> temp(s);
        while(true) {
            temp /= rank;
            if(temp < BasicConst::getInstance().getExpectedRelativeError())
                break;
            result += temp;
            temp *= s;
            ++rank;
        }
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> exp(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedScalarUnit>(1));
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(1));
        Scalar<MultiPrecision, true> temp(s);
        while(true) {
            temp /= rank;
            if(temp < BasicConst::getInstance().getExpectedRelativeError())
                break;
            result += temp;
            temp *= s;
            ++rank;
        }
        return result;
    }

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
    Scalar<MultiPrecision, false> cos(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedScalarUnit>(1));
        if(s == BasicConst::getInstance().get_0())
            return result;
        Scalar<MultiPrecision, false> square_n = square(s);
        Scalar<MultiPrecision, false> temp_1(square_n);
        Scalar<MultiPrecision, false> temp_2(static_cast<SignedScalarUnit>(2));
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(2));
        bool changeSign = true;

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= s;
            rank += BasicConst::getInstance().get_1();
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_n;
            temp_2 *= rank;
            rank += BasicConst::getInstance().get_1();
            temp_2 *= rank;
        }
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> cos(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedScalarUnit>(1));
        if(s == BasicConst::getInstance().get_0())
            return result;
        Scalar<MultiPrecision, true> square_n = square(s);
        Scalar<MultiPrecision, true> temp_1(square_n);
        Scalar<MultiPrecision, true> temp_2(static_cast<SignedScalarUnit>(2));
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(2));
        bool changeSign = true;

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= s;
            rank += BasicConst::getInstance().get_1();
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_n;
            temp_2 *= rank;
            rank += BasicConst::getInstance().get_1();
            temp_2 *= rank;
        }
        return result;
    }

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
    Scalar<MultiPrecision, false> sin(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedScalarUnit>(0));
        if(s == BasicConst::getInstance().get_0())
            return result;
        Scalar<MultiPrecision, false> square_s = square(s);
        Scalar<MultiPrecision, false> temp_1(s);
        Scalar<MultiPrecision, false> temp_2(static_cast<SignedScalarUnit>(1));
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(1));
        bool changeSign = false;

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= s;
            rank += BasicConst::getInstance().get_1();
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_s;
            temp_2 *= rank;
            rank += BasicConst::getInstance().get_1();
            temp_2 *= rank;
        }
        return Scalar<MultiPrecision, false>(result);
    }

    template<>
    Scalar<MultiPrecision, true> sin(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedScalarUnit>(0));
        if(s == BasicConst::getInstance().get_0())
            return result;
        Scalar<MultiPrecision, true> square_s = square(s);
        Scalar<MultiPrecision, true> temp_1(s);
        Scalar<MultiPrecision, true> temp_2(static_cast<SignedScalarUnit>(1));
        Scalar<MultiPrecision, false> rank(static_cast<SignedScalarUnit>(1));
        bool changeSign = false;

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= s;
            rank += BasicConst::getInstance().get_1();
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_s;
            temp_2 *= rank;
            rank += BasicConst::getInstance().get_1();
            temp_2 *= rank;
        }
        return Scalar<MultiPrecision, true>(result);
    }

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
    Scalar<MultiPrecision, false> arccos(const Scalar<MultiPrecision, false>& s) {
        return Solve::bisectionMethod(cos, s, BasicConst::getInstance().get_0(), MathConst::getInstance().getPI()
                , BasicConst::getInstance().get_1(), BasicConst::getInstance().getMinus_1());
    }

    template<>
    Scalar<MultiPrecision, true> arccos(const Scalar<MultiPrecision, true>& s) {
        return Solve::bisectionMethod(cos, s, BasicConst::getInstance().get_0(), MathConst::getInstance().getPI()
                , BasicConst::getInstance().get_1(), BasicConst::getInstance().getMinus_1());
    }

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
    Scalar<MultiPrecision, false> arcsin(const Scalar<MultiPrecision, false>& s) {
        return Solve::bisectionMethod(sin, s, MathConst::getInstance().getMinus_PI_2(), MathConst::getInstance().getPI_2()
                , BasicConst::getInstance().getMinus_1(), BasicConst::getInstance().get_1());
    }

    template<>
    Scalar<MultiPrecision, true> arcsin(const Scalar<MultiPrecision, true>& s) {
        return Solve::bisectionMethod(sin, s, MathConst::getInstance().getMinus_PI_2(), MathConst::getInstance().getPI_2()
                , BasicConst::getInstance().getMinus_1(), BasicConst::getInstance().get_1());
    }

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
    Scalar<MultiPrecision, true> arctan(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result = arcsin(s / sqrt(square(s) + BasicConst::getInstance().get_1()));
        if((result.getLength() ^ s.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

    template<>
    Scalar<MultiPrecision, false> arctan(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result = arcsin(s / sqrt(square(s) + BasicConst::getInstance().get_1()));
        if((result.getLength() ^ s.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

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
    Scalar<MultiPrecision, false> cosh(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result = exp(s);
        result = (result + reciprocal(result)) >> 1;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> cosh(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result = exp(s);
        result = (result + reciprocal(result)) >> 1;
        return result;
    }

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
    Scalar<MultiPrecision, false> sinh(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result = exp(s);
        Scalar<MultiPrecision, false> temp = reciprocal(result);
        result -= temp;
        result >>= 1;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> sinh(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result = exp(s);
        Scalar<MultiPrecision, true> temp = reciprocal(result);
        result -= temp;
        result >>= 1;
        return result;
    }

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
    Scalar<MultiPrecision, false> tanh(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result = exp(s);
        Scalar<MultiPrecision, false> temp = reciprocal(result);
        Scalar<MultiPrecision, false> temp1 = result + temp;
        result -= temp;
        result /= temp1;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> tanh(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result = exp(s);
        Scalar<MultiPrecision, true> temp = reciprocal(result);
        Scalar<MultiPrecision, true> temp1 = result + temp;
        result -= temp;
        result /= temp1;
        return result;
    }

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
    Scalar<MultiPrecision, false> sech(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedScalarUnit>(2));
        Scalar<MultiPrecision, false> temp = exp(s);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> sech(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedScalarUnit>(2));
        Scalar<MultiPrecision, true> temp = exp(s);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

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
    Scalar<MultiPrecision, false> csch(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedScalarUnit>(2));
        Scalar<MultiPrecision, false> temp = exp(s);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> csch(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedScalarUnit>(2));
        Scalar<MultiPrecision, true> temp = exp(s);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

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
    Scalar<MultiPrecision, false> coth(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result = exp(s);
        Scalar<MultiPrecision, false> temp = reciprocal(result);
        Scalar<MultiPrecision, false> temp1 = result - temp;
        result += temp;
        result /= temp1;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> coth(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result = exp(s);
        Scalar<MultiPrecision, true> temp = reciprocal(result);
        Scalar<MultiPrecision, true> temp1 = result - temp;
        result += temp;
        result /= temp1;
        return result;
    }

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
    Scalar<MultiPrecision, false> arccosh(const Scalar<MultiPrecision, false>& s) {
        return ln(sqrt(square(s) - BasicConst::getInstance().get_1()) + s);
    }

    template<>
    Scalar<MultiPrecision, true> arccosh(const Scalar<MultiPrecision, true>& s) {
        return ln(sqrt(square(s) - BasicConst::getInstance().get_1()) + s);
    }

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
    Scalar<MultiPrecision, false> arcsinh(const Scalar<MultiPrecision, false>& s) {
        return ln(sqrt(square(s) + BasicConst::getInstance().get_1()) + s);
    }

    template<>
    Scalar<MultiPrecision, true> arcsinh(const Scalar<MultiPrecision, true>& s) {
        return ln(sqrt(square(s) + BasicConst::getInstance().get_1()) + s);
    }

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
    Scalar<MultiPrecision, false> arctanh(const Scalar<MultiPrecision, false>& s) {
        return ln((BasicConst::getInstance().get_1() + s) / (BasicConst::getInstance().get_1() - s)) >> 1;
    }

    template<>
    Scalar<MultiPrecision, true> arctanh(const Scalar<MultiPrecision, true>& s) {
        return ln((BasicConst::getInstance().get_1() + s) / (BasicConst::getInstance().get_1() - s)) >> 1;
    }

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
    Scalar<MultiPrecision, false> arccoth(const Scalar<MultiPrecision, false>& s) {
        return ln((s + BasicConst::getInstance().get_1()) / (s - BasicConst::getInstance().get_1())) >> 1;
    }

    template<>
    Scalar<MultiPrecision, true> arccoth(const Scalar<MultiPrecision, true>& s) {
        return ln((s + BasicConst::getInstance().get_1()) / (s - BasicConst::getInstance().get_1())) >> 1;
    }
}

#endif