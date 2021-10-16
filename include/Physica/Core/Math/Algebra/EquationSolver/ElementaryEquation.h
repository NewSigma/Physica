/*
 * Copyright 2020-2021 WeiBo He.
 *
 * This file is part of Physica.

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
#pragma once

namespace Physica::Core {
    /**
     * Equations are defined like this:
     * ScalarType func(const ScalarType&)
     */
    template<class Function, class ScalarType>
    ScalarType bisectionMethod(
            Function func,
            const ScalarType& x1,
            const ScalarType& x2) {
        const ScalarType y1 = func(x1);
        const ScalarType y2 = func(x2);
        return bisectionMethod<Function, ScalarType>(func, x1, x2, y1, y2);
    }

    template<class Function, class ScalarType>
    ScalarType bisectionMethod(
            Function func,
            const ScalarType& x1,
            const ScalarType& x2,
            const ScalarType& y1,
            const ScalarType& y2) {
        if(y1.isZero())
            return x1;
        if(y2.isZero())
            return x2;
        assert(!ScalarType::matchSign(y1, y2)); //Root must be existent

        ScalarType result = (x1 + x2) >> 1;
        ScalarType y_result = ScalarType::One();

        ScalarType error = ScalarType(x1 - x2).toAbs() >> 1;
        ScalarType x_left(x1);
        ScalarType x_right(x2);
        ScalarType y_left(y1);

        bool delta_left_sign = !y_left.isPositive();
        bool delta_right_sign;
        do {
            y_result = func(result);
            delta_right_sign = !y_result.isPositive();

            if(delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign =!y_left.isPositive();
            }
            else
                x_right = result;
            result = (x_left + x_right) >> 1;
            error >>= 1;
        } while(result.getPower() - error.getPower() < GlobalPrecision);
        result.toUnitA();
        return result;
    }
}
