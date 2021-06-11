/*
 * Copyright 2021 WeiBo He.
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
#pragma once

#include <cassert>

namespace Physica::Core {
    template<class T, class Function>
    T bisection(Function func, const T& n, const T& x1, const T& x2, const T& y1, const T& y2) {
        constexpr double epsilon_trivial = std::numeric_limits<typename T::TrivialType>::epsilon();
        const T epsilon = T(epsilon_trivial);
        assert(T(n - y1).isPositive() ^ T(n - y2).isPositive()); //(n - y1) and (n - y2) have different sign
        if(n == y1)
            return T(x1);
        if(n == y2)
            return T(x2);

        T result = (x1 + x2) / T(2);

        T error = (x1 - x2) / T(2);
        error.toAbs();
        T x_left(x1);
        T x_right(x2);
        T y_left(y1);

        T y_result;
        bool delta_left_sign = n > y_left;
        bool delta_right_sign;
        do {
            y_result = func(result);
            delta_right_sign = n > y_result;

            if(delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign = n > y_left;
            }
            else
                x_right = result;
            result = (x_left + x_right) / T(2);
            error /= T(2);
        } while(error > abs(epsilon * result));
        return T(result).toUnitA();
    }

    template<class T, class Function>
    T bisection(Function func, const T& n, const T& x1, const T& x2) {
        T y1 = func(x1);
        T y2 = func(x2);
        return bisection(func, n, x1, x2, y1, y2);
    }
}
