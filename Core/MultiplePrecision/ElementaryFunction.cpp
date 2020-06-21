/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <climits>
#include "Core/Header/Scalar.h"
#include "Core/Header/Solve.h"

namespace Physica::Core {
    //Return a real number between 0 and 1.
    Scalar randomNumerical() {
        return Scalar(double(random()) / RAND_MAX);
    }
    //Return a real number lowerBound and upperBound.
    Scalar randomNumerical(const Scalar& lowerBound, const Scalar& upperBound) {
        return randomNumerical() * (upperBound - lowerBound) + lowerBound;
    }
    //Reference: GMP Doc BaseCase Multiplication
    Scalar square(const Scalar& n) {
        if(n == BasicConst::getInstance().get_1())
            return Scalar(n);
        else {
            auto n_size = n.getSize();
            //Estimate the ed of result first. we will calculate it accurately later.
            const auto length = 2 * n_size;
            Scalar result(length, n.power * 2 + 1);

            for(int i = 0; i < n_size - 1; ++i)
                result[i + n_size] = mulAddArrByWord(result.byte + i + i + 1
                        , n.byte + i + 1, n_size - i - 1, n.byte[i]);
            //Shift count is known, possible to optimize the performance.
            byteLeftShiftEq(result.byte, length, 1);

            ScalarUnit high, low;
            for(int i = 0; i < n_size; ++i) {
                mulWordByWord(high, low, n.byte[i], n.byte[i]);
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

    Scalar floor(const Scalar& n) {
        if(n.isInteger())
            return Scalar(n);
        const auto size = n.getSize();
        const auto power = n.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = n.isNegative() ? -length : length;
        Scalar result(length, power);
        for(int i = 0; i < length; ++i)
            result[i] = n[i];
        return result;
    }

    Scalar ceil(const Scalar& n) {
        auto f = floor(n);
        return ++f;
    }

    Scalar reciprocal(const Scalar& n) {
        return BasicConst::getInstance().get_1() / n;
    }
    /*
     * *_light functions do not consider the error caused by a. For example, sqrt_light does not calculate
     * (sqrt(n + a) - sqrt(n)) for error.
     */
    Scalar sqrt_light(const Scalar& n) {
        if(n.isNegative())
            qFatal("Can not resolve the square root of a minus value.");
        if(n.isZero())
            return Scalar(BasicConst::getInstance().get_0());
        Scalar copy_n(n);
        //Let n < 1 so as to control error.
        int add_power = 0;
        if(copy_n.getPower() > 0) {
            if(copy_n.getPower() % 2 == 0) {
                add_power = copy_n.getPower() / 2 + 1;
                copy_n.power = -2;
            }
            else {
                add_power = (copy_n.getPower() + 1) / 2;
                copy_n.power = -1;
            }
        }

        Scalar result = getOne();
        //3.33 is the big approximate value of ln(10)/ln(2)
        for(int i = 0; i < LONG_WIDTH * BasicConst::getInstance().GlobalPrecision; ++i)
            result = (result + Scalar::div(copy_n, result)) >> 1U;
        result.power += add_power;
        result.toUnitA();

        return result;
    }

    Scalar sqrt(const Scalar& n) {
        Scalar result  = sqrt_light(n);
        if(n.getA() != 0) {
            Scalar n_error = n.getMinimum();
            if(n_error.isNegative())
                n_error = n.getMaximum();
            Scalar error = sqrt_light(n_error);
            error -= result;
            error.toAbs();

            result.applyError(error);
        }
        return result;
    }
//TODO not completed: Use gamma function.
    Scalar factorial(const Scalar& n) {
        if(n.getLength() < 0)
            qFatal("Can not resolve the factorial of a minus value.");
        if(!n.isInteger())
            qFatal("Can not resolve the factorial of a float value.");

        Scalar result = getOne();
        Scalar temp = getOne();
        while(temp < n) {
            temp += BasicConst::getInstance().get_1();
            result *= temp;
        }
        return result;
    }

    Scalar ln_light(const Scalar& n) {
        if(!n.isPositive())
            qFatal("Can not resolve the logarithm of zero or a negative value.");
        if(n == BasicConst::getInstance().get_1())
            return getZero();
        Scalar result = getZero();
        Scalar temp_1 = Scalar::sub(n, BasicConst::getInstance().get_1())
                        / Scalar::add(n, BasicConst::getInstance().get_1());
        Scalar copy_temp_1(temp_1);
        Scalar rank = getOne();

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
            if(result.getPower() - criteria.getPower() >= BasicConst::getInstance().GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= copy_temp_1;
            rank += BasicConst::getInstance().get_1();
        }
        result *= BasicConst::getInstance().get_2();
        return result;
    }

    Scalar ln(const Scalar& n) {
        Scalar result = ln_light(n);
        if(n.getA() != 0) {
            Scalar n_error = n.getMinimum();
            if(n_error.isNegative())
                n_error = n.getMaximum();
            Scalar error = ln_light(n_error);
            error -= result;
            error.toAbs();

            result.applyError(error);
        }
        return result;
    }
    //Return log_a n
    Scalar log(const Scalar& n, const Scalar& a) {
        if(!n.isPositive() || !a.isPositive())
            qFatal("Can not resolve the logarithm of zero or a negative value.");
        return ln(n) / ln(a);
    }

    Scalar exp(const Scalar& n) {
        Scalar result = getOne();
        Scalar temp(n);
        Scalar rank = getOne();
        while(true) {
            temp /= rank;
            if(temp < BasicConst::getInstance().getExpectedRelativeError())
                break;
            result += temp;
            temp *= n;
            ++rank;
        }
        return result;
    }
    /*
     * Taylor's formula n.th term: (-1)^n * x^(2n) / (2n!)
     * Here temp_1 = x^(2n), temp_2 = 2n!, rank = 2n
     */
    Scalar cos(const Scalar& n) {
        Scalar result = getOne();
        if(n == BasicConst::getInstance().get_0())
            return result;
        Scalar square_n = square(n);
        Scalar temp_1(square_n);
        Scalar temp_2 = getTwo();
        Scalar rank = getTwo();
        bool changeSign = true;

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= n;
            rank += BasicConst::getInstance().get_1();
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= BasicConst::getInstance().GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_n;
            temp_2 *= rank;
            rank += BasicConst::getInstance().get_1();
            temp_2 *= rank;
        }
        return result;
    }

    Scalar sin(const Scalar& n) {
        Scalar result = getZero();
        if(n == BasicConst::getInstance().get_0())
            return result;
        Scalar square_n = square(n);
        Scalar temp_1(n);
        Scalar temp_2 = getOne();
        Scalar rank = getOne();
        bool changeSign = false;

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= n;
            rank += BasicConst::getInstance().get_1();
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= BasicConst::getInstance().GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_n;
            temp_2 *= rank;
            rank += BasicConst::getInstance().get_1();
            temp_2 *= rank;
        }
        return result;
    }

    Scalar tan(const Scalar& n) {
        return sin(n) / cos(n);
    }

    Scalar sec(const Scalar& n) {
        return reciprocal(cos(n));
    }

    Scalar csc(const Scalar& n) {
        return reciprocal(sin(n));
    }

    Scalar cot(const Scalar& n) {
        return cos(n) / sin(n);
    }

    Scalar arccos(const Scalar& n) {
        return Solve::bisectionMethod(cos, n, BasicConst::getInstance().get_0(), MathConst::getInstance().getPI()
                , BasicConst::getInstance().get_1(), BasicConst::getInstance().getMinus_1());
    }

    Scalar arcsin(const Scalar& n) {
        return Solve::bisectionMethod(sin, n, MathConst::getInstance().getMinus_PI_2(), MathConst::getInstance().getPI_2()
                , BasicConst::getInstance().getMinus_1(), BasicConst::getInstance().get_1());
    }

    Scalar arctan(const Scalar& n) {
        Scalar result = arcsin(n / sqrt(square(n) + BasicConst::getInstance().get_1()));
        if((result.getLength() ^ n.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

    Scalar arcsec(const Scalar& n) {
        return arccos(reciprocal(n));
    }

    Scalar arccsc(const Scalar& n) {
        return arcsin(reciprocal(n));
    }

    Scalar arccot(const Scalar& n) {
        return arctan(reciprocal(n));
    }

    Scalar cosh(const Scalar& n) {
        Scalar result = exp(n);
        result = (result + reciprocal(result)) / BasicConst::getInstance().get_2();
        return result;
    }

    Scalar sinh(const Scalar& n) {
        Scalar result = exp(n);
        Scalar temp = reciprocal(result);
        result -= temp;
        result /= BasicConst::getInstance().get_2();
        return result;
    }

    Scalar tanh(const Scalar& n) {
        Scalar result = exp(n);
        Scalar temp = reciprocal(result);
        Scalar temp1 = result + temp;
        result -= temp;
        result /= temp1;
        return result;
    }

    Scalar sech(const Scalar& n) {
        Scalar result = getTwo();
        Scalar temp = exp(n);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    Scalar csch(const Scalar& n) {
        Scalar result = getTwo();
        Scalar temp = exp(n);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

    Scalar coth(const Scalar& n) {
        Scalar result = exp(n);
        Scalar temp = reciprocal(result);
        Scalar temp1 = result - temp;
        result += temp;
        result /= temp1;
        return result;
    }

    Scalar arccosh(const Scalar& n) {
        return ln(sqrt(square(n) - BasicConst::getInstance().get_1()) + n);
    }

    Scalar arcsinh(const Scalar& n) {
        return ln(sqrt(square(n) + BasicConst::getInstance().get_1()) + n);
    }

    Scalar arctanh(const Scalar& n) {
        return ln((BasicConst::getInstance().get_1() + n) / (BasicConst::getInstance().get_1() - n)) / BasicConst::getInstance().get_2();
    }

    Scalar arcsech(const Scalar& n) {
        return arccosh(reciprocal(n));
    }

    Scalar arccsch(const Scalar& n) {
        return arcsinh(reciprocal(n));
    }

    Scalar arccoth(const Scalar& n) {
        return ln((n + BasicConst::getInstance().get_1()) / (n - BasicConst::getInstance().get_1())) / BasicConst::getInstance().get_2();
    }
}