/*
 * Copyright 2020-2023 Weibo He.
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
    Scalar<MultiPrecision> abs(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> temp(s);
        return Scalar<MultiPrecision>(std::move(temp.toAbs()));
    }

    template<>
    Scalar<MultiPrecision> square(const Scalar<MultiPrecision>& s) noexcept {
        if(s == BasicConst::getInstance()._1)
            return Scalar(s);
        else {
            const auto s_size = s.getSize();
            //Estimate the length of result first. we will calculate it accurately later.
            auto length = 2 * s_size;
            Scalar<MultiPrecision> result(length, s.getPower() * 2 + 1);

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

    template<>
    Scalar<MultiPrecision> sqrt(const Scalar<MultiPrecision>& s) noexcept {
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

        Scalar result = Scalar<MultiPrecision>(static_cast<SignedMPUnit>(1));
        for(unsigned int i = 0; i < MPUnitWidth * GlobalPrecision; ++i)
            result = (result + Scalar<MultiPrecision>::div(copy_s, result)) >> 1U;
        result.power += add_power;
        return result;
    }

    template<>
    Scalar<MultiPrecision> ln(const Scalar<MultiPrecision>& s) noexcept {
        assert(s.isPositive());
        Scalar<MultiPrecision> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._1)
            return result;
        const auto& _1 = BasicConst::getInstance()._1;
        auto temp_1 = Scalar<MultiPrecision>::sub(s, _1)
                      / Scalar<MultiPrecision>::add(s, _1);
        Scalar<MultiPrecision> copy_temp_1(temp_1);
        Scalar<MultiPrecision> rank(static_cast<SignedMPUnit>(1));

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
    Scalar<MultiPrecision> exp(const Scalar<MultiPrecision>& s) noexcept {
        if(s.isNegative())
            return reciprocal(exp(-s));
        Scalar<MultiPrecision> result = 1;
        Scalar<MultiPrecision> rank = 1;
        Scalar<MultiPrecision> temp(s);
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
    Scalar<MultiPrecision> cos(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result(static_cast<SignedMPUnit>(1));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<MultiPrecision> square_n = square(s);
        Scalar<MultiPrecision> temp_1(square_n);
        Scalar<MultiPrecision> temp_2(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision> rank(static_cast<SignedMPUnit>(2));
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
    Scalar<MultiPrecision> sin(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result(static_cast<SignedMPUnit>(0));
        if(s == BasicConst::getInstance()._0)
            return result;
        Scalar<MultiPrecision> square_s = square(s);
        Scalar<MultiPrecision> temp_1(s);
        Scalar<MultiPrecision> temp_2(static_cast<SignedMPUnit>(1));
        Scalar<MultiPrecision> rank(static_cast<SignedMPUnit>(1));
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
        return Scalar<MultiPrecision>(result);
    }

    template<>
    Scalar<MultiPrecision> arccos(const Scalar<MultiPrecision>& s) noexcept {
        using ScalarType = Scalar<MultiPrecision>;
        auto func = [&](const ScalarType& x) -> ScalarType { return cos(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           BasicConst::getInstance()._0,
                                                           MathConst::getInstance().PI,
                                                           BasicConst::getInstance()._1,
                                                           BasicConst::getInstance().Minus_1);
    }

    template<>
    Scalar<MultiPrecision> arcsin(const Scalar<MultiPrecision>& s) noexcept {
        using ScalarType = Scalar<MultiPrecision>;
        auto func = [&](const ScalarType& x) -> ScalarType { return sin(x) - s; };
        return bisectionMethod<decltype(func), ScalarType>(func,
                                                           MathConst::getInstance().Minus_PI_2,
                                                           MathConst::getInstance().PI_2,
                                                           BasicConst::getInstance().Minus_1,
                                                           BasicConst::getInstance()._1);
    }

    template<>
    Scalar<MultiPrecision> arctan(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> temp = square(s) + BasicConst::getInstance()._1;
        Scalar<MultiPrecision> result = arcsin(s / sqrt(temp));
        if((result.getLength() ^ s.getLength()) < 0) // NOLINT(hicpp-signed-bitwise)
            result.toAbs();
        return result;
    }

    template<>
    Scalar<MultiPrecision> cosh(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result = exp(s);
        result = (result + reciprocal(result)) >> 1;
        return result;
    }

    template<>
    Scalar<MultiPrecision> sinh(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result = exp(s);
        Scalar<MultiPrecision> temp = reciprocal(result);
        result -= temp;
        result >>= 1;
        return result;
    }

    template<>
    Scalar<MultiPrecision> tanh(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result = exp(s);
        Scalar<MultiPrecision> temp = reciprocal(result);
        Scalar<MultiPrecision> temp1 = result + temp;
        result -= temp;
        result /= temp1;
        return result;
    }

    template<>
    Scalar<MultiPrecision> sech(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision> temp = exp(s);
        temp += reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<MultiPrecision> csch(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result(static_cast<SignedMPUnit>(2));
        Scalar<MultiPrecision> temp = exp(s);
        temp -= reciprocal(temp);
        result /= temp;
        return result;
    }

    template<>
    Scalar<MultiPrecision> coth(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> result = exp(s);
        Scalar<MultiPrecision> temp = reciprocal(result);
        Scalar<MultiPrecision> temp1 = result - temp;
        result += temp;
        result /= temp1;
        return result;
    }

    template<>
    Scalar<MultiPrecision> arccosh(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> temp = square(s) - BasicConst::getInstance()._1;
        Scalar<MultiPrecision> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<MultiPrecision> arcsinh(const Scalar<MultiPrecision>& s) noexcept {
        Scalar<MultiPrecision> temp = square(s) + BasicConst::getInstance()._1;
        Scalar<MultiPrecision> temp1 = sqrt(temp) + s;
        return ln(temp1);
    }

    template<>
    Scalar<MultiPrecision> arctanh(const Scalar<MultiPrecision>& s) noexcept {
        return ln((BasicConst::getInstance()._1 + s)
                  / Scalar<MultiPrecision>(BasicConst::getInstance()._1 - s)) >> 1;
    }

    template<>
    Scalar<MultiPrecision> arccoth(const Scalar<MultiPrecision>& s) noexcept {
        return ln((s + BasicConst::getInstance()._1)
                  / Scalar<MultiPrecision>(s - BasicConst::getInstance()._1)) >> 1;
    }
}