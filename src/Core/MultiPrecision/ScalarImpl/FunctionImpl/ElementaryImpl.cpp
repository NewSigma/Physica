/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/MultiPrecition/ScalarImpl/FunctionImpl/BisectionMethod.h"

namespace Physica::Core {
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
    /*!
     * Implemented using taylor's series: e^x = 1 + x + x^2/2 + x^3/3! + ... + x^n/n! + e^(ax)*x^(n+1)/(n+1)!  (0 < a < 1)
     * We use the first n terms and have a error: E = e^(ax)*x^(n+1)/(n+1)!
     *
     * If x > 0, we have E = e^(ax)*x^(n+1)/(n+1)! < e^(x)*x^(n+1)/(n+1)!
     * So the relative error is E/e^x < x^(n+1)/(n+1)!
     * If x < 0, we calculate reciprocal of e^-x instead.
     */
    template<>
    Scalar<MultiPrecision, false> exp(const Scalar<MultiPrecision, false>& s) {
        if(s.isNegative())
            return reciprocal(exp(-s));
        Scalar<MultiPrecision, false> result = 1;
        Scalar<MultiPrecision, false> rank = 1;
        Scalar<MultiPrecision, false> temp(s);
        const auto& relativeError = BasicConst::getInstance().getExpectedRelativeError();
        while(true) {
            temp /= rank;
            if(absCompare(relativeError, temp))
                break;
            result += temp;
            temp *= s;
            ++rank;
        }
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> exp(const Scalar<MultiPrecision, true>& s) {
        if(s.isNegative())
            return reciprocal(exp(-s));
        Scalar<MultiPrecision, true> result = 1;
        Scalar<MultiPrecision, false> rank = 1;
        Scalar<MultiPrecision, true> temp(s);
        const auto& relativeError = BasicConst::getInstance().getExpectedRelativeError();
        int a = 0;
        while(true) {
            ++a;
            temp /= rank;
            if(absCompare(relativeError, temp))
                break;
            result += temp;
            temp *= s;
            ++rank;
        }
        return result;
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

    template<>
    Scalar<MultiPrecision, false> arccos(const Scalar<MultiPrecision, false>& s) {
        return bisectionMethod<MultiPrecision, false>(cos, s
                , BasicConst::getInstance().get_0(), MathConst::getInstance().getPI()
                , BasicConst::getInstance().get_1(), BasicConst::getInstance().getMinus_1());
    }

    template<>
    Scalar<MultiPrecision, true> arccos(const Scalar<MultiPrecision, true>& s) {
        return bisectionMethod<MultiPrecision, true>(cos, s
                , BasicConst::getInstance().get_0(), MathConst::getInstance().getPI()
                , BasicConst::getInstance().get_1(), BasicConst::getInstance().getMinus_1());
    }

    template<>
    Scalar<MultiPrecision, false> arcsin(const Scalar<MultiPrecision, false>& s) {
        return bisectionMethod<MultiPrecision, false>(sin, s
                , MathConst::getInstance().getMinus_PI_2(), MathConst::getInstance().getPI_2()
                , BasicConst::getInstance().getMinus_1(), BasicConst::getInstance().get_1());
    }

    template<>
    Scalar<MultiPrecision, true> arcsin(const Scalar<MultiPrecision, true>& s) {
        return bisectionMethod<MultiPrecision, true>(sin, s
                , MathConst::getInstance().getMinus_PI_2(), MathConst::getInstance().getPI_2()
                , BasicConst::getInstance().getMinus_1(), BasicConst::getInstance().get_1());
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

    template<>
    Scalar<MultiPrecision, false> arccosh(const Scalar<MultiPrecision, false>& s) {
        return ln(sqrt(square(s) - BasicConst::getInstance().get_1()) + s);
    }

    template<>
    Scalar<MultiPrecision, true> arccosh(const Scalar<MultiPrecision, true>& s) {
        return ln(sqrt(square(s) - BasicConst::getInstance().get_1()) + s);
    }

    template<>
    Scalar<MultiPrecision, false> arcsinh(const Scalar<MultiPrecision, false>& s) {
        return ln(sqrt(square(s) + BasicConst::getInstance().get_1()) + s);
    }

    template<>
    Scalar<MultiPrecision, true> arcsinh(const Scalar<MultiPrecision, true>& s) {
        return ln(sqrt(square(s) + BasicConst::getInstance().get_1()) + s);
    }

    template<>
    Scalar<MultiPrecision, false> arctanh(const Scalar<MultiPrecision, false>& s) {
        return ln((BasicConst::getInstance().get_1() + s) / (BasicConst::getInstance().get_1() - s)) >> 1;
    }

    template<>
    Scalar<MultiPrecision, true> arctanh(const Scalar<MultiPrecision, true>& s) {
        return ln((BasicConst::getInstance().get_1() + s) / (BasicConst::getInstance().get_1() - s)) >> 1;
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