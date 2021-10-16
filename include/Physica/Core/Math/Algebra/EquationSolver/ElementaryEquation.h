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

        const ScalarType half = ScalarType(0.5);
        ScalarType result = (x1 + x2) * half;
        ScalarType y_result = ScalarType::One();

        ScalarType error = ScalarType(x1 - x2).toAbs() * half;
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
            result = (x_left + x_right) * half;
            error *= half;
        } while(abs(result * std::numeric_limits<ScalarType>::epsilon()) < error);
        result.toUnitA();
        return result;
    }
    /**
     * Reference:
     * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:559-560
     */
    template<class Function, class ScalarType>
    ScalarType secant(
            Function func,
            const ScalarType& x1,
            const ScalarType& x2,
            const ScalarType& abs_error) {
        const ScalarType y1 = func(x1);
        const ScalarType y2 = func(x2);
        return secant<Function, ScalarType>(func, x1, x2, y1, y2, abs_error);
    }

    template<class Function, class ScalarType>
    ScalarType secant(
            Function func,
            const ScalarType& x1,
            const ScalarType& x2,
            const ScalarType& y1,
            const ScalarType& y2,
            const ScalarType& abs_error) {
        if(y1.isZero())
            return x1;
        if(y2.isZero())
            return x2;
        assert(!ScalarType::matchSign(y1, y2)); //Root must be existent

        ScalarType x_old = x1;
        ScalarType x_now = x2;
        ScalarType y_old = y1;
        ScalarType y_now = y2;
        do {
            ScalarType temp = (x_old * y_now - x_now * y_old) / (y_now - y_old);
            x_old = std::move(x_now);
            y_old = std::move(y_now);
            x_now = std::move(temp);
            y_now = func(x_now);
        } while(abs(x_now - x_old) > abs_error);
        return x_now;
    }
}
