/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <climits>
#include "Core/Header/Numerical.h"
#include "Core/Header/Solve.h"

namespace Physica::Core {
    //Return a real number between 0 and 1.
    Numerical randomNumerical() {
        return Numerical(double(random()) / RAND_MAX);
    }
    //Return a real number lowerBound and upperBound.
    Numerical randomNumerical(const Numerical& lowerBound, const Numerical& upperBound) {
        return randomNumerical() * (upperBound - lowerBound) + lowerBound;
    }
    //Reference: GMP Doc BaseCase Multiplication
    Numerical square(const Numerical& n) {
        if(n == BasicConst::getInstance().get_1())
            return Numerical(n);
        else {
            auto n_size = n.getSize();
            //Estimate the ed of result first. we will calculate it accurately later.
            const auto length = 2 * n_size;
            Numerical result(length, n.power * 2 + 1);

            for(int i = 0; i < n_size - 1; ++i)
                result[i + n_size] = mulAddArrByWord(result.byte + i + i + 1
                        , n.byte + i + 1, n_size - i - 1, n.byte[i]);
            //Shift count is known, possible to optimize the performance.
            byteLeftShiftEq(result.byte, length, 1);

            NumericalUnit high, low;
            for(int i = 0; i < n_size; ++i) {
                mulWordByWord(high, low, n.byte[i], n.byte[i]);
                result[2 * i] += low;
                result[2 * i + 1] += high;
            }

            if(high == 0) {
                --result.power;
                result.byte =
                        reinterpret_cast<NumericalUnit*>(realloc(result.byte, (length - 1) * sizeof(NumericalUnit)));
            }
            return result;
        }
    }

    Numerical floor(const Numerical& n) {
        if(n.isInteger())
            return Numerical(n);
        const auto size = n.getSize();
        const auto power = n.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = n.isNegative() ? -length : length;
        Numerical result(length, power);
        for(int i = 0; i < length; ++i)
            result[i] = n[i];
        return result;
    }

    Numerical ceil(const Numerical& n) {
        auto f = floor(n);
        return ++f;
    }

    Numerical reciprocal(const Numerical& n) {
        return BasicConst::getInstance().get_1() / n;
    }
    /*
     * *_light functions do not consider the error caused by a. For example, sqrt_light does not calculate
     * (sqrt(n + a) - sqrt(n)) for error.
     */
    Numerical sqrt_light(const Numerical& n) {
        if(n.isNegative())
            qFatal("Can not resolve the square root of a minus value.");
        if(n.isZero())
            return Numerical(BasicConst::getInstance().get_0());
        Numerical copy_n(n);
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

        Numerical result = getOne();
        //3.33 is the big approximate value of ln(10)/ln(2)
        for(int i = 0; i < LONG_WIDTH * BasicConst::getInstance().GlobalPrecision; ++i)
            result = (result + Numerical::div(copy_n, result)) >> 1U;
        result.power += add_power;
        result.toUnitA();

        return result;
    }

    Numerical sqrt(const Numerical& n) {
        Numerical result  = sqrt_light(n);
        if(n.getA() != 0) {
            Numerical n_error = n.getMinimum();
            if(n_error.isNegative())
                n_error = n.getMaximum();
            Numerical error = sqrt_light(n_error);
            error -= result;
            error.toAbs();

            result.applyError(error);
        }
        return result;
    }
//TODO not completed: Use gamma function.
    Numerical factorial(const Numerical& n) {
        if(n.getLength() < 0)
            qFatal("Can not resolve the factorial of a minus value.");
        if(!n.isInteger())
            qFatal("Can not resolve the factorial of a float value.");

        Numerical result = getOne();
        Numerical temp = getOne();
        while(temp < n) {
            temp += BasicConst::getInstance().get_1();
            result *= temp;
        }
        return result;
    }

    Numerical ln_light(const Numerical& n) {
        if(!n.isPositive())
            qFatal("Can not resolve the logarithm of zero or a negative value.");
        if(n == BasicConst::getInstance().get_1())
            return getZero();
        Numerical result = getZero();
        Numerical temp_1 = Numerical::sub(n, BasicConst::getInstance().get_1())
                           / Numerical::add(n, BasicConst::getInstance().get_1());
        Numerical copy_temp_1(temp_1);
        Numerical rank = getOne();

        temp_1.toUnitA();
        copy_temp_1.toUnitA();
        while(true) {
            //Calculate one term of the taylor series.
            Numerical temp = temp_1 / rank;
            temp.clearA();
            result += temp;

            temp_1 *= copy_temp_1;
            rank += BasicConst::getInstance().get_1();
            Numerical criteria = temp_1 / rank;
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

    Numerical ln(const Numerical& n) {
        Numerical result = ln_light(n);
        if(n.getA() != 0) {
            Numerical n_error = n.getMinimum();
            if(n_error.isNegative())
                n_error = n.getMaximum();
            Numerical error = ln_light(n_error);
            error -= result;
            error.toAbs();

            result.applyError(error);
        }
        return result;
    }
    //Return log_a n
    Numerical log(const Numerical& n, const Numerical& a) {
        if(!n.isPositive() || !a.isPositive())
            qFatal("Can not resolve the logarithm of zero or a negative value.");
        return ln(n) / ln(a);
    }

    Numerical exp(const Numerical& n) {
        Numerical result = getOne();
        Numerical temp(n);
        Numerical rank = getOne();
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
    Numerical cos(const Numerical& n) {
        Numerical result = getOne();
        if(n == BasicConst::getInstance().get_0())
            return result;
        Numerical square_n = square(n);
        Numerical temp_1(square_n);
        Numerical temp_2 = getTwo();
        Numerical rank = getTwo();
        bool changeSign = true;

        while(true) {
            //Calculate one term of the taylor series.
            Numerical temp = temp_1 / temp_2;
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

    Numerical sin(const Numerical& n) {
        Numerical result = getZero();
        if(n == BasicConst::getInstance().get_0())
            return result;
        Numerical square_n = square(n);
        Numerical temp_1(n);
        Numerical temp_2 = getOne();
        Numerical rank = getOne();
        bool changeSign = false;

        while(true) {
            //Calculate one term of the taylor series.
            Numerical temp = temp_1 / temp_2;
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

    Numerical tan(const Numerical& n) {
        return sin(n) / cos(n);
    }

    Numerical sec(const Numerical& n) {
        return reciprocal(cos(n));
    }

    Numerical csc(const Numerical& n) {
        return reciprocal(sin(n));
    }

    Numerical cot(const Numerical& n) {
        return cos(n) / sin(n);
    }

    Numerical arccos(const Numerical& n) {
        return Solve::bisectionMethod(cos, n, BasicConst::getInstance().get_0(), MathConst::getInstance().getPI()
                , BasicConst::getInstance().get_1(), BasicConst::getInstance().getMinus_1());
    }

    Numerical arcsin(const Numerical& n) {
        return Solve::bisectionMethod(sin, n, MathConst::getInstance().getMinus_PI_2(), MathConst::getInstance().getPI_2()
                , BasicConst::getInstance().getMinus_1(), BasicConst::getInstance().get_1());
    }

    Numerical arctan(const Numerical& n) {
        Numerical result = arcsin(n / sqrt(square(n) + BasicConst::getInstance().get_1()));
        if((result.getLength() ^ n.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

    Numerical arcsec(const Numerical& n) {
        return arccos(reciprocal(n));
    }

    Numerical arccsc(const Numerical& n) {
        return arcsin(reciprocal(n));
    }

    Numerical arccot(const Numerical& n) {
        return arctan(reciprocal(n));
    }

    Numerical cosh(const Numerical& n) {
        Numerical result = exp(n);
        result = (result + reciprocal(result)) / BasicConst::getInstance().get_2();
        return result;
    }

    Numerical sinh(const Numerical& n) {
        Numerical result = exp(n);
        Numerical temp = reciprocal(result);
        result -= temp;
        result /= BasicConst::getInstance().get_2();
        return result;
    }

    Numerical tanh(const Numerical& n) {
        Numerical result = exp(n);
        Numerical temp = reciprocal(result);
        Numerical temp1 = result + temp;
        result -= temp;
        result /= temp1;
        return result;
    }

    Numerical sech(const Numerical& n) {
        Numerical result = getTwo();
        Numerical temp = exp(n);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    Numerical csch(const Numerical& n) {
        Numerical result = getTwo();
        Numerical temp = exp(n);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

    Numerical coth(const Numerical& n) {
        Numerical result = exp(n);
        Numerical temp = reciprocal(result);
        Numerical temp1 = result - temp;
        result += temp;
        result /= temp1;
        return result;
    }

    Numerical arccosh(const Numerical& n) {
        return ln(sqrt(square(n) - BasicConst::getInstance().get_1()) + n);
    }

    Numerical arcsinh(const Numerical& n) {
        return ln(sqrt(square(n) + BasicConst::getInstance().get_1()) + n);
    }

    Numerical arctanh(const Numerical& n) {
        return ln((BasicConst::getInstance().get_1() + n) / (BasicConst::getInstance().get_1() - n)) / BasicConst::getInstance().get_2();
    }

    Numerical arcsech(const Numerical& n) {
        return arccosh(reciprocal(n));
    }

    Numerical arccsch(const Numerical& n) {
        return arcsinh(reciprocal(n));
    }

    Numerical arccoth(const Numerical& n) {
        return ln((n + BasicConst::getInstance().get_1()) / (n - BasicConst::getInstance().get_1())) / BasicConst::getInstance().get_2();
    }
}