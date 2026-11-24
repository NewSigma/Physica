/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Scalar/Real.h"

namespace Physica {
    template<Scalar T>
    T bisectionMethod(std::invocable<T> auto fn, const T& x1, const T& x2) {
        const T y1 = fn(x1);
        const T y2 = fn(x2);
        return bisectionMethod<T>(fn, x1, x2, y1, y2);
    }

    template<Scalar T>
    T bisectionMethod(std::invocable<T> auto fn, const T& x1, const T& x2, const T& y1, const T& y2) {
        if (y1.isZero())
            return x1;
        if (y2.isZero())
            return x2;
        assert(!T::matchSign(y1, y2)); // Root must be existent

        const T half = T(0.5);
        T result = (x1 + x2) * half;
        T y_result(1);

        T error = abs(x1 - x2) * half;
        T x_left(x1);
        T x_right(x2);
        T y_left(y1);

        bool delta_left_sign = !y_left.isPositive();
        bool delta_right_sign;
        do {
            y_result = fn(result);
            delta_right_sign = !y_result.isPositive();

            if (delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign = !y_left.isPositive();
            }
            else
                x_right = result;
            result = (x_left + x_right) * half;
            error *= half;
        } while (abs(result * std::numeric_limits<T>::epsilon()) < error);
        return result;
    }
    /**
     * Reference:
     * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:559-560
     */
    template<Scalar T>
    T secant(std::invocable<T> auto fn, const T& x1, const T& x2, const T& abs_error) {
        const T y1 = fn(x1);
        const T y2 = fn(x2);
        return secant<T>(fn, x1, x2, y1, y2, abs_error);
    }

    template<Scalar T>
    T secant(std::invocable<T> auto fn, const T& x1, const T& x2, const T& y1, const T& y2, const T& abs_error) {
        if (y1.isZero())
            return x1;
        if (y2.isZero())
            return x2;
        assert(!T::matchSign(y1, y2)); // Root must be existent

        T x_old = x1;
        T x_now = x2;
        T y_old = y1;
        T y_now = y2;
        do {
            T temp = (x_old * y_now - x_now * y_old) / (y_now - y_old);
            x_old = std::move(x_now);
            y_old = std::move(y_now);
            x_now = std::move(temp);
            y_now = fn(x_now);
        } while (abs(x_now - x_old) > abs_error);
        return x_now;
    }
}
