/*
 * Copyright 2020-2026 Weibo He.
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
#include "Physica/Core/Math/Calculus/ElementaryEquation.h"

namespace Physica {
    template<>
    Real<FloatMP> abs(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> temp(s);
        return Real<FloatMP>(std::move(temp.toAbs()));
    }

    template<>
    Real<FloatMP> sqrt(const Real<FloatMP>& s) noexcept {
        assert(!s.isNegative());
        if(s.isZero())
            return Real(BasicConst::getInstance()._0);
        Real copy_s(s);
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

        Real result = Real<FloatMP>(static_cast<SignedMPUnit>(1));
        for(unsigned int i = 0; i < MPUnitWidth * Real<FloatMP>::GlobalPrecision; ++i)
            result = (result + Real<FloatMP>::div(copy_s, result)) >> 1U;
        result.power += add_power;
        return result;
    }

    template<>
    Real<FloatMP> ln(const Real<FloatMP>& s) noexcept {
        assert(s.isPositive());
        Real<FloatMP> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._1)
            return result;
        const auto& _1 = BasicConst::getInstance()._1;
        auto temp_1 = Real<FloatMP>::sub(s, _1)
                      / Real<FloatMP>::add(s, _1);
        Real<FloatMP> copy_temp_1(temp_1);
        Real<FloatMP> rank(static_cast<SignedMPUnit>(1));

        while(true) {
            //Calculate one term of the taylor series.
            Real temp = temp_1 / rank;
            result += temp;

            temp_1 *= copy_temp_1;
            rank += BasicConst::getInstance()._1;
            Real criteria = temp_1 / rank;
            //Break if result meets the precision goal.
            if(result.getPower() - criteria.getPower() >= Real<FloatMP>::GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= copy_temp_1;
            rank += _1;
        }
        result *= BasicConst::getInstance()._2;
        return result;
    }

    template<>
    Real<FloatMP> exp(const Real<FloatMP>& s) noexcept {
        if(s.isNegative())
            return reciprocal(exp(-s));
        Real<FloatMP> result = 1;
        Real<FloatMP> rank = 1;
        Real<FloatMP> temp(s);
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
    Real<FloatMP> pow(const Real<FloatMP>& s1, const Real<FloatMP>& s2) noexcept;

    template<>
    Real<FloatMP> factorial(const Real<FloatMP>& s) noexcept {
        //Optimize: Unnecessary copy during floor() if s is a integer itself.
        const Real<FloatMP> integer = floor(s);

        Real<FloatMP> result(SignedMPUnit(1));
        Real<FloatMP> temp(SignedMPUnit(1));
        while(temp < integer)
            result *= ++temp;
        return result;
    }

    template<>
    Real<FloatMP> cos(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result(static_cast<SignedMPUnit>(1));
        if(s == BasicConst::getInstance()._0)
            return result;
        Real<FloatMP> square_n = square(s);
        Real<FloatMP> temp_1(square_n);
        Real<FloatMP> temp_2(static_cast<SignedMPUnit>(2));
        Real<FloatMP> rank(static_cast<SignedMPUnit>(2));
        bool changeSign = true;

        while(true) {
            //Calculate one term of the taylor series.
            Real temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= s;
            rank += BasicConst::getInstance()._1;
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= Real<FloatMP>::GlobalPrecision)
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
    Real<FloatMP> sin(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._0)
            return result;
        Real<FloatMP> square_s = square(s);
        Real<FloatMP> temp_1(s);
        Real<FloatMP> temp_2(static_cast<SignedMPUnit>(1));
        Real<FloatMP> rank(static_cast<SignedMPUnit>(1));
        bool changeSign = false;

        while(true) {
            //Calculate one term of the taylor series.
            Real temp = temp_1 / temp_2;
            if(changeSign)
                temp.toOpposite();
            changeSign = !changeSign;
            result += temp;
            //Here the temp means the criteria of break.
            temp *= s;
            rank += BasicConst::getInstance()._1;
            temp /= rank;
            //Break if result meets the precision goal.
            if(result.getPower() - temp.getPower() >= Real<FloatMP>::GlobalPrecision)
                break;
            //Prepare for next calculate.
            temp_1 *= square_s;
            temp_2 *= rank;
            rank += BasicConst::getInstance()._1;
            temp_2 *= rank;
        }
        return Real<FloatMP>(result);
    }

    template<>
    Real<FloatMP> arccos(const Real<FloatMP>& s) noexcept {
        using T = Real<FloatMP>;
        auto func = [&](const T& x) -> T { return cos(x) - s; };
        return bisectionMethod<T>(
                func,
                BasicConst::getInstance()._0,
                MathConst<T>::getInstance().PI,
                BasicConst::getInstance()._1,
                BasicConst::getInstance().Minus_1);
    }

    template<>
    Real<FloatMP> arcsin(const Real<FloatMP>& s) noexcept {
        using T = Real<FloatMP>;
        auto func = [&](const T& x) -> T { return sin(x) - s; };
        return bisectionMethod<T>(
                func,
                MathConst<T>::getInstance().Minus_PI_2,
                MathConst<T>::getInstance().PI_2,
                BasicConst::getInstance().Minus_1,
                BasicConst::getInstance()._1);
    }

    template<>
    Real<FloatMP> arctan(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> temp = square(s) + BasicConst::getInstance()._1;
        Real<FloatMP> result = arcsin(s / sqrt(temp));
        if((result.getLength() ^ s.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

    template<>
    Real<FloatMP> cosh(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result = exp(s);
        result = (result + reciprocal(result)) >> 1;
        return result;
    }

    template<>
    Real<FloatMP> sinh(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result = exp(s);
        Real<FloatMP> temp = reciprocal(result);
        result -= temp;
        result >>= 1;
        return result;
    }

    template<>
    Real<FloatMP> tanh(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result = exp(s);
        Real<FloatMP> temp = reciprocal(result);
        Real<FloatMP> temp1 = result + temp;
        result -= temp;
        result /= temp1;
        return result;
    }

    template<>
    Real<FloatMP> sech(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result(static_cast<SignedMPUnit>(2));
        Real<FloatMP> temp = exp(s);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Real<FloatMP> csch(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result(static_cast<SignedMPUnit>(2));
        Real<FloatMP> temp = exp(s);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Real<FloatMP> coth(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> result = exp(s);
        Real<FloatMP> temp = reciprocal(result);
        Real<FloatMP> temp1 = result - temp;
        result += temp;
        result /= temp1;
        return result;
    }

    template<>
    Real<FloatMP> arccosh(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> temp = square(s) - BasicConst::getInstance()._1;
        Real<FloatMP> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Real<FloatMP> arcsinh(const Real<FloatMP>& s) noexcept {
        Real<FloatMP> temp = square(s) + BasicConst::getInstance()._1;
        Real<FloatMP> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Real<FloatMP> arctanh(const Real<FloatMP>& s) noexcept {
        return ln((BasicConst::getInstance()._1 + s)
                  / Real<FloatMP>(BasicConst::getInstance()._1 - s)) >> 1;
    }

    template<>
    Real<FloatMP> arccoth(const Real<FloatMP>& s) noexcept {
        return ln((s + BasicConst::getInstance()._1)
                  / Real<FloatMP>(s - BasicConst::getInstance()._1)) >> 1;
    }

    template<>
    Real<FloatMP> floor(const Real<FloatMP>& s) noexcept {
        if(s.isInteger())
            return Real(s);
        const auto size = s.getSize();
        const auto power = s.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = s.isNegative() ? -length : length;
        Real<FloatMP> result(length, power);
        for(int i = 0; i < length; ++i)
            result.setByte(i, s[i]);
        return result;
    }
}
