/*
 * Copyright 2021-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:142
     */
    template<Scalar T>
    void chebyshev_fit(const T& from, const T& to, Vector auto& coeff, std::invocable<T> auto fn) {
        assert(from < to);
        const size_t n = coeff.getLength();
        const T rep = reciprocal(T(n));
        const auto funcArr = VectorND<T>::generate([from, to, fn, rep](size_t i) {
            const T temp1 = (to - from) * T(0.5);
            const T temp2 = (to + from) * T(0.5);
            const T y = cospi((T(i) + 0.5) * rep);
            return fn(fma(y, temp1, temp2));
        }, n);

        const T factor = rep * 2;
        for (size_t i = 0; i < n; ++i) {
            T sum(0);
            for (size_t j = 0; j < n; ++j)
                sum += funcArr[j] * cospi(T(i) * (T(j) + 0.5) * rep);
            coeff[i] = factor * sum;
        }
    }
    /**
     * Fit a even function, the performance is better than chebyshev_fit()
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:143
     */
    template<Scalar T>
    void chebyshev_fit_even(const T& from, const T& to, Vector auto& coeff, std::invocable<T> auto fn) {
        assert(from < to);
        const size_t n = coeff.getLength();
        const size_t n2 = n * 2;
        const T rep = reciprocal(T(n2));
        const auto funcArr = VectorND<T>::generate([from, to, fn, rep](size_t i) {
            const T temp1 = (to - from) * T(0.5);
            const T temp2 = (to + from) * T(0.5);
            const T y = cospi((T(i) + 0.5) * rep);
            return fn(fma(y, temp1, temp2));
        }, n2);

        const T factor = rep * 4;
        for (size_t i = 0; i < n; ++i) {
            T sum(0);
            for (size_t j = 0; j < n; ++j)
                sum += funcArr[j] * cospi(T(2 * i) * (T(j) + 0.5) * rep);
            coeff[i] = factor * sum;
        }
    }
    /**
     * Fit a odd function, the performance is better than chebyshev_fit()
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:143
     */
    template<Scalar T>
    void chebyshev_fit_odd(const T& from, const T& to, Vector auto& coeff, std::invocable<T> auto fn) {
        assert(from < to);
        const size_t n = coeff.getLength();
        const size_t n2 = n * 2;
        const T rep = reciprocal(T(n2));
        const auto funcArr = VectorND<T>::generate([from, to, fn, rep](size_t i) {
            const T temp1 = (to - from) * T(0.5);
            const T temp2 = (to + from) * T(0.5);
            const T y = cospi((T(i) + 0.5) * rep);
            const T x = fma(y, temp1, temp2);
            return fn(x) / x;
        }, n2);

        const T factor = rep * 4;
        for (size_t i = 0; i < n; ++i) {
            T sum(0);
            for (size_t j = 0; j < n; ++j)
                sum += funcArr[j] * cospi(T(2 * i) * (T(j) + 0.5) * rep);
            coeff[i] = factor * sum;
        }
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:142
     */
    template<Scalar T>
    T chebyshev_calc(const T& from, const T& to, const Vector auto& coeff, const T& x) {
        assert(from <= x && x <= to);
        const T y = T(fma(x, T(2), -(to + from))) / T(to - from);
        const T y2 = y * T(2);
        T d1(0), d2(0);
        for (size_t i = coeff.getLength() - 1; i > 0; --i) {
            const T temp = d1;
            d1 = fma(y2, d1, coeff.calc(i) - d2);
            d2 = temp;
        }
        return fma(y, d1, fma(coeff.calc(0), T(0.5), -d2));
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:143
     */
    template<Scalar T>
    T chebyshev_calc_even(const T& from, const T& to, const Vector auto& coeff, const T& x) {
        assert(from <= x && x <= to);
        return chebyshev_calc<T>(from, to, coeff, fma(x * 2, x, T(-1)));
    }

    template<Scalar T>
    T chebyshev_calc_odd(const T& from, const T& to, const Vector auto& coeff, const T& x) {
        assert(from <= x && x <= to);
        return chebyshev_calc_even<T>(from, to, coeff, x) * x;
    }
}
