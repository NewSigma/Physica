/*
 * Copyright 2020-2021 WeiBo He.
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
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/EquationSolver/ElementaryEquation.h"

namespace Physica::Core {
    template<>
    Scalar<MultiPrecision, false> sqrt(const Scalar<MultiPrecision, false>& s) {
        assert(!s.isNegative());
        if(s.isZero())
            return Scalar(BasicConst::getInstance()._0);
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

        Scalar result = Scalar<MultiPrecision, false>(static_cast<SignedMPUnit>(1));
        for(unsigned int i = 0; i < MPUnitWidth * GlobalPrecision; ++i)
            result = (result + Scalar<MultiPrecision, false>::divNoError(copy_s, result)) >> 1U;
        result.power += add_power;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> sqrt(const Scalar<MultiPrecision, true>& s) {
        assert(!s.isNegative());
        if(s.isZero())
            return Scalar<MultiPrecision, true>(BasicConst::getInstance()._0);
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

        Scalar result = Scalar<MultiPrecision, true>(static_cast<SignedMPUnit>(1));
        for(unsigned int i = 0; i < MPUnitWidth * GlobalPrecision; ++i)
            result = (result + Scalar<MultiPrecision, false>::divNoError(copy_s, result)) >> 1U;
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
        assert(s.isPositive());
        Scalar<MultiPrecision, false> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._1)
            return result;
        const auto& _1 = BasicConst::getInstance()._1;
        auto temp_1 = Scalar<MultiPrecision, false>::subNoError(s, _1)
                      / Scalar<MultiPrecision, false>::addNoError(s, _1);
        Scalar<MultiPrecision, false> copy_temp_1(temp_1);
        Scalar<MultiPrecision, false> rank(static_cast<SignedMPUnit>(1));

        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / rank;
            result += temp;

            temp_1 *= copy_temp_1;
            rank += BasicConst::getInstance()._1;
            Scalar criteria = temp_1 / rank;
            //Break if result meets the precision goal.
            if(result.getPower() - criteria.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= copy_temp_1;
            rank += _1;
        }
        result *= BasicConst::getInstance()._2;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> ln(const Scalar<MultiPrecision, true>& s) {
        assert(s.isPositive());
        Scalar<MultiPrecision, true> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._1)
            return result;
        const auto& _1 = BasicConst::getInstance()._1;
        auto temp_1 = Scalar<MultiPrecision, false>::subWithError(s, _1)
                      / Scalar<MultiPrecision, false>::addNoError(s, _1);
        Scalar<MultiPrecision, true> copy_temp_1(temp_1);
        Scalar<MultiPrecision, false> rank(static_cast<SignedMPUnit>(1));

        temp_1.toUnitA();
        copy_temp_1.toUnitA();
        while(true) {
            //Calculate one term of the taylor series.
            Scalar temp = temp_1 / rank;
            temp.clearA();
            result += temp;

            temp_1 *= copy_temp_1;
            rank += BasicConst::getInstance()._1;
            Scalar criteria = temp_1 / rank;
            //Break if result meets the precision goal.
            if(result.getPower() - criteria.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= copy_temp_1;
            rank += _1;
        }
        result *= BasicConst::getInstance()._2;

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
        const auto& relativeError = BasicConst::getInstance().expectedRelativeError;
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
        Scalar<MultiPrecision, true> result(1);
        Scalar<MultiPrecision, false> rank(1);
        Scalar<MultiPrecision, true> temp(s);
        const auto& relativeError = BasicConst::getInstance().expectedRelativeError;
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
    Scalar<MultiPrecision, false> cos(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedMPUnit>(1));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<MultiPrecision, false> square_n = square(s);
        Scalar<MultiPrecision, false> temp_1(square_n);
        Scalar<MultiPrecision, false> temp_2(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision, false> rank(static_cast<SignedMPUnit>(2));
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
            rank += BasicConst::getInstance()._1;
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_n;
            temp_2 *= rank;
            rank += BasicConst::getInstance()._1;
            temp_2 *= rank;
        }
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> cos(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedMPUnit>(1));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<MultiPrecision, true> square_n = square(s);
        Scalar<MultiPrecision, true> temp_1(square_n);
        Scalar<MultiPrecision, true> temp_2(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision, false> rank(static_cast<SignedMPUnit>(2));
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
            rank += BasicConst::getInstance()._1;
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_n;
            temp_2 *= rank;
            rank += BasicConst::getInstance()._1;
            temp_2 *= rank;
        }
        return result;
    }

    template<>
    Scalar<MultiPrecision, false> sin(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<MultiPrecision, false> square_s = square(s);
        Scalar<MultiPrecision, false> temp_1(s);
        Scalar<MultiPrecision, false> temp_2(static_cast<SignedMPUnit>(1));
        Scalar<MultiPrecision, false> rank(static_cast<SignedMPUnit>(1));
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
            rank += BasicConst::getInstance()._1;
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_s;
            temp_2 *= rank;
            rank += BasicConst::getInstance()._1;
            temp_2 *= rank;
        }
        return Scalar<MultiPrecision, false>(result);
    }

    template<>
    Scalar<MultiPrecision, true> sin(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<MultiPrecision, true> square_s = square(s);
        Scalar<MultiPrecision, true> temp_1(s);
        Scalar<MultiPrecision, true> temp_2(static_cast<SignedMPUnit>(1));
        Scalar<MultiPrecision, false> rank(static_cast<SignedMPUnit>(1));
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
            rank += BasicConst::getInstance()._1;
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_s;
            temp_2 *= rank;
            rank += BasicConst::getInstance()._1;
            temp_2 *= rank;
        }
        return Scalar<MultiPrecision, true>(result);
    }

    template<>
    Scalar<MultiPrecision, false> arccos(const Scalar<MultiPrecision, false>& s) {
        using ScalarType = Scalar<MultiPrecision, false>;
        auto func = [&](const ScalarType& x) -> ScalarType { return cos(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           BasicConst::getInstance()._0,
                                                           MathConst::getInstance().PI,
                                                           BasicConst::getInstance()._1,
                                                           BasicConst::getInstance().Minus_1);
    }

    template<>
    Scalar<MultiPrecision, true> arccos(const Scalar<MultiPrecision, true>& s) {
        using ScalarType = Scalar<MultiPrecision, true>;
        auto func = [&](const ScalarType& x) -> ScalarType { return cos(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           BasicConst::getInstance()._0,
                                                           MathConst::getInstance().PI,
                                                           BasicConst::getInstance()._1,
                                                           BasicConst::getInstance().Minus_1);
    }

    template<>
    Scalar<MultiPrecision, false> arcsin(const Scalar<MultiPrecision, false>& s) {
        using ScalarType = Scalar<MultiPrecision, false>;
        auto func = [&](const ScalarType& x) -> ScalarType { return sin(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           MathConst::getInstance().Minus_PI_2,
                                                           MathConst::getInstance().PI_2,
                                                           BasicConst::getInstance().Minus_1,
                                                           BasicConst::getInstance()._1);
    }

    template<>
    Scalar<MultiPrecision, true> arcsin(const Scalar<MultiPrecision, true>& s) {
        using ScalarType = Scalar<MultiPrecision, true>;
        auto func = [&](const ScalarType& x) -> ScalarType { return sin(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           MathConst::getInstance().Minus_PI_2,
                                                           MathConst::getInstance().PI_2,
                                                           BasicConst::getInstance().Minus_1,
                                                           BasicConst::getInstance()._1);
    }
    //!FixIt: arctan() does not consider accuracy.
    template<>
    Scalar<MultiPrecision, true> arctan(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> temp = square(s) + BasicConst::getInstance()._1;
        Scalar<MultiPrecision, true> result = arcsin(s / sqrt(temp));
        if((result.getLength() ^ s.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

    template<>
    Scalar<MultiPrecision, false> arctan(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> temp = square(s) + BasicConst::getInstance()._1;
        Scalar<MultiPrecision, false> result = arcsin(s / sqrt(temp));
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
        Scalar<MultiPrecision, false> result(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision, false> temp = exp(s);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> sech(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision, true> temp = exp(s);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<MultiPrecision, false> csch(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> result(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision, false> temp = exp(s);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<MultiPrecision, true> csch(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> result(static_cast<SignedMPUnit>(2));
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
        Scalar<MultiPrecision, false> temp = square(s) - BasicConst::getInstance()._1;
        Scalar<MultiPrecision, false> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<MultiPrecision, true> arccosh(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> temp = square(s) - BasicConst::getInstance()._1;
        Scalar<MultiPrecision, true> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<MultiPrecision, false> arcsinh(const Scalar<MultiPrecision, false>& s) {
        Scalar<MultiPrecision, false> temp = square(s) + BasicConst::getInstance()._1;
        Scalar<MultiPrecision, false> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<MultiPrecision, true> arcsinh(const Scalar<MultiPrecision, true>& s) {
        Scalar<MultiPrecision, true> temp = square(s) - BasicConst::getInstance()._1;
        Scalar<MultiPrecision, true> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<MultiPrecision, false> arctanh(const Scalar<MultiPrecision, false>& s) {
        return ln((BasicConst::getInstance()._1 + s)
                  / Scalar<MultiPrecision, false>(BasicConst::getInstance()._1 - s)) >> 1;
    }

    template<>
    Scalar<MultiPrecision, true> arctanh(const Scalar<MultiPrecision, true>& s) {
        return ln((BasicConst::getInstance()._1 + s)
                  / Scalar<MultiPrecision, true>(BasicConst::getInstance()._1 - s)) >> 1;
    }

    template<>
    Scalar<MultiPrecision, false> arccoth(const Scalar<MultiPrecision, false>& s) {
        return ln((s + BasicConst::getInstance()._1)
                  / Scalar<MultiPrecision, false>(s - BasicConst::getInstance()._1)) >> 1;
    }

    template<>
    Scalar<MultiPrecision, true> arccoth(const Scalar<MultiPrecision, true>& s) {
        return ln((s + BasicConst::getInstance()._1)
                  / Scalar<MultiPrecision, true>(s - BasicConst::getInstance()._1)) >> 1;
    }
}