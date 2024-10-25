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
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/EquationSolver/ElementaryEquation.h"

namespace Physica::Core {
    template<>
    Scalar<FloatMP> abs(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> temp(s);
        return Scalar<FloatMP>(std::move(temp.toAbs()));
    }

    template<>
    Scalar<FloatMP> sqrt(const Scalar<FloatMP>& s) noexcept {
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

        Scalar result = Scalar<FloatMP>(static_cast<SignedMPUnit>(1));
        for(unsigned int i = 0; i < MPUnitWidth * GlobalPrecision; ++i)
            result = (result + Scalar<FloatMP>::div(copy_s, result)) >> 1U;
        result.power += add_power;
        return result;
    }

    template<>
    Scalar<FloatMP> ln(const Scalar<FloatMP>& s) noexcept {
        assert(s.isPositive());
        Scalar<FloatMP> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._1)
            return result;
        const auto& _1 = BasicConst::getInstance()._1;
        auto temp_1 = Scalar<FloatMP>::sub(s, _1)
                      / Scalar<FloatMP>::add(s, _1);
        Scalar<FloatMP> copy_temp_1(temp_1);
        Scalar<FloatMP> rank(static_cast<SignedMPUnit>(1));

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
    Scalar<FloatMP> exp(const Scalar<FloatMP>& s) noexcept {
        if(s.isNegative())
            return reciprocal(exp(-s));
        Scalar<FloatMP> result = 1;
        Scalar<FloatMP> rank = 1;
        Scalar<FloatMP> temp(s);
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
    Scalar<FloatMP> pow(const Scalar<FloatMP>& s1, const Scalar<FloatMP>& s2) noexcept;

    template<>
    Scalar<FloatMP> factorial(const Scalar<FloatMP>& s) noexcept {
        //Optimize: Unnecessary copy during floor() if s is a integer itself.
        const Scalar<FloatMP> integer = floor(s);

        Scalar<FloatMP> result(SignedMPUnit(1));
        Scalar<FloatMP> temp(SignedMPUnit(1));
        while(temp < integer)
            result *= ++temp;
        return result;
    }

    template<>
    Scalar<FloatMP> cos(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result(static_cast<SignedMPUnit>(1));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<FloatMP> square_n = square(s);
        Scalar<FloatMP> temp_1(square_n);
        Scalar<FloatMP> temp_2(static_cast<SignedMPUnit>(2));
        Scalar<FloatMP> rank(static_cast<SignedMPUnit>(2));
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
    Scalar<FloatMP> sin(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<FloatMP> square_s = square(s);
        Scalar<FloatMP> temp_1(s);
        Scalar<FloatMP> temp_2(static_cast<SignedMPUnit>(1));
        Scalar<FloatMP> rank(static_cast<SignedMPUnit>(1));
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
        return Scalar<FloatMP>(result);
    }

    template<>
    Scalar<FloatMP> arccos(const Scalar<FloatMP>& s) noexcept {
        using ScalarType = Scalar<FloatMP>;
        auto func = [&](const ScalarType& x) -> ScalarType { return cos(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           BasicConst::getInstance()._0,
                                                           MathConst::getInstance().PI,
                                                           BasicConst::getInstance()._1,
                                                           BasicConst::getInstance().Minus_1);
    }

    template<>
    Scalar<FloatMP> arcsin(const Scalar<FloatMP>& s) noexcept {
        using ScalarType = Scalar<FloatMP>;
        auto func = [&](const ScalarType& x) -> ScalarType { return sin(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           MathConst::getInstance().Minus_PI_2,
                                                           MathConst::getInstance().PI_2,
                                                           BasicConst::getInstance().Minus_1,
                                                           BasicConst::getInstance()._1);
    }

    template<>
    Scalar<FloatMP> arctan(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> temp = square(s) + BasicConst::getInstance()._1;
        Scalar<FloatMP> result = arcsin(s / sqrt(temp));
        if((result.getLength() ^ s.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

    template<>
    Scalar<FloatMP> cosh(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result = exp(s);
        result = (result + reciprocal(result)) >> 1;
        return result;
    }

    template<>
    Scalar<FloatMP> sinh(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result = exp(s);
        Scalar<FloatMP> temp = reciprocal(result);
        result -= temp;
        result >>= 1;
        return result;
    }

    template<>
    Scalar<FloatMP> tanh(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result = exp(s);
        Scalar<FloatMP> temp = reciprocal(result);
        Scalar<FloatMP> temp1 = result + temp;
        result -= temp;
        result /= temp1;
        return result;
    }

    template<>
    Scalar<FloatMP> sech(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result(static_cast<SignedMPUnit>(2));
        Scalar<FloatMP> temp = exp(s);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<FloatMP> csch(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result(static_cast<SignedMPUnit>(2));
        Scalar<FloatMP> temp = exp(s);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<FloatMP> coth(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> result = exp(s);
        Scalar<FloatMP> temp = reciprocal(result);
        Scalar<FloatMP> temp1 = result - temp;
        result += temp;
        result /= temp1;
        return result;
    }

    template<>
    Scalar<FloatMP> arccosh(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> temp = square(s) - BasicConst::getInstance()._1;
        Scalar<FloatMP> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<FloatMP> arcsinh(const Scalar<FloatMP>& s) noexcept {
        Scalar<FloatMP> temp = square(s) + BasicConst::getInstance()._1;
        Scalar<FloatMP> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<FloatMP> arctanh(const Scalar<FloatMP>& s) noexcept {
        return ln((BasicConst::getInstance()._1 + s)
                  / Scalar<FloatMP>(BasicConst::getInstance()._1 - s)) >> 1;
    }

    template<>
    Scalar<FloatMP> arccoth(const Scalar<FloatMP>& s) noexcept {
        return ln((s + BasicConst::getInstance()._1)
                  / Scalar<FloatMP>(s - BasicConst::getInstance()._1)) >> 1;
    }

    template<>
    Scalar<FloatMP> floor(const Scalar<FloatMP>& s) noexcept {
        if(s.isInteger())
            return Scalar(s);
        const auto size = s.getSize();
        const auto power = s.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = s.isNegative() ? -length : length;
        Scalar<FloatMP> result(length, power);
        for(int i = 0; i < length; ++i)
            result.setByte(i, s[i]);
        return result;
    }
}