/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Math/MathConst.h"

namespace Physica::Core {
    /*!
     * Consts that need some calculates.
     * Should call new to Const_1 so as to make calculates available.
     */
    MathConst<FloatMP>::MathConst() {
        //0.31 is the big approximation of ln(2) / ln(10)
        PI = T(calcPI(
                static_cast<int>(static_cast<double>(MPUnitWidth) * GlobalPrecision * 0.31) + 1));
        E = T(exp(BasicConst::getInstance()._1));

        PI_2 = T(PI >> 1);
        Minus_PI_2 = T(-PI_2);
    }
    /*!
     * precision is the number of effective digits in decimal.
     * Reference:
     * [1] http://www.pi314.net/eng/salamin.php
     * [2] https://blog.csdn.net/liangbch/article/details/78724041
     */
    MathConst<FloatMP>::T MathConst<FloatMP>::calcPI(int precision) {
        const auto& basicConst = BasicConst::getInstance();

        T a(static_cast<SignedMPUnit>(1));
        T x(static_cast<SignedMPUnit>(1));
        T b(reciprocal(sqrt(basicConst._2)));
        T c(reciprocal(basicConst._4));

        int goal = 1;
        while(goal < precision) {
            Real y(a);
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

    const MathConst<FloatMP>& MathConst<FloatMP>::getInstance() noexcept {
        static MathConst mathConst{};
        return mathConst;
    }
}
