/*
 * Copyright 2019 WeiBo He.
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

namespace Physica::Core {
    /*!
     * Basic consts that initialize directly.
     */
    BasicConst::BasicConst()
            : ln_2(std::log(2))
            , ln_10(std::log(10))
            , ln_2_10(std::log(2) / std::log(10))
            , plotPoints(static_cast<SignedMPUnit>(20)) {
        expectedRelativeError = Scalar<MultiPrecision, false>(1, 1 - GlobalPrecision);
        expectedRelativeError.setByte(0, 1);

        stepSize = Scalar<MultiPrecision, false>(1, - GlobalPrecision / 2); //(- GlobalPrecision / 2) is selected by experience.
        stepSize.setByte(0, 1);

        R_MAX = 2147483647;
        _0 = 0;
        _1 = 1;
        Minus_1 = -1;
        _2 = 2;
        Minus_2 = -2;
        _3 = 3;
        Minus_3 = -3;
        _4 = 4;
        Minus_4 = -4;
        _10 = 10;
    }
    /*!
     * Consts that need some calculates.
     * Should call new to Const_1 so as to make calculates available.
     */
    MathConst::MathConst() {
        //0.31 is the big approximation of ln(2) / ln(10)
        PI = MultiScalar(calcPI(
                static_cast<int>(static_cast<double>(MPUnitWidth) * GlobalPrecision * 0.31) + 1));
        E = MultiScalar(exp(BasicConst::getInstance()._1));

        PI_2 = MultiScalar(PI >> 1);
        Minus_PI_2 = MultiScalar(-PI_2);
    }
    /*!
     * precision is the number of effective digits in decimal.
     * Reference:
     * http://www.pi314.net/eng/salamin.php
     * https://blog.csdn.net/liangbch/article/details/78724041
     */
    MultiScalar MathConst::calcPI(int precision) {
        const auto& basicConst = BasicConst::getInstance();

        MultiScalar a(static_cast<SignedMPUnit>(1));
        MultiScalar x(static_cast<SignedMPUnit>(1));
        MultiScalar b(reciprocal(sqrt(basicConst._2)));
        MultiScalar c(reciprocal(basicConst._4));

        int goal = 1;
        while(goal < precision) {
            Scalar y(a);
            a = (a + b) >> 1;
            b = sqrt(b * y);
            y -= a;
            c -= y * y * x;
            x *= basicConst._2;
            goal *= 2;
        }
        a = (a + b) >> 1;
        return a * a / c;
    }
}
