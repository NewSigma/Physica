/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    template<Scalar T, class Functor>
    [[nodiscard]] T bisection(Functor func, const T& target, const T& x1, const T& x2, const T& y1, const T& y2) noexcept(std::is_nothrow_invocable_v<Functor, T>) {
        const T epsilon = std::numeric_limits<T>::epsilon();
        assert(T(target - y1).isPositive() ^ T(target - y2).isPositive() && "[Error]: (target - y1) and (target - y2) must have different sign");
        if (target == y1)
            return T(x1);
        if (target == y2)
            return T(x2);

        T result = (x1 + x2) / T(2);

        T error = abs((x1 - x2) * T(0.5));
        T x_left(x1);
        T x_right(x2);
        T y_left(y1);

        T y_result;
        bool delta_left_sign = target > y_left;
        bool delta_right_sign;
        do {
            y_result = func(result);
            delta_right_sign = target > y_result;

            if (delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign = target > y_left;
            }
            else
                x_right = result;
            result = (x_left + x_right) / T(2);
            error /= T(2);
        } while (error > abs(epsilon * result));
        return result;
    }

    template<Scalar T, class Functor>
    [[nodiscard]] T bisection(Functor func, const T& target, const T& x1, const T& x2) noexcept(std::is_nothrow_invocable_v<Functor, T>) {
        return bisection(func, target, x1, x2, func(x1), func(x2));
    }
}
