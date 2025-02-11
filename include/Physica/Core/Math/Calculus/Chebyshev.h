/*
 * Copyright 2021-2024 Weibo He.
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

#include <cmath>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:142
     */
    template<Scalar T, Vector U, class Function>
    void chebyshev_fit(const T& from, const T& to, LValueVector<U>& coeff, Function func) {
        constexpr double pi = M_PI;
        assert(from < to);

        const size_t n = coeff.getLength();
        const double n_1 = 1.0 / n;
        VectorND<T> funcArr(n);
        /* Fill func arr */ {
            const T& temp1 = (to - from) * T(0.5);
            const T& temp2 = (to + from) * T(0.5);
            for (size_t i = 0; i < n; ++i) {
                const T y = T(std::cos(pi * (i + 0.5) * n_1)); //Optimize: 1. pi * n_1 can be stored in a variable  2. use add to avoid mul  3. store pi / 2 and calculate 2 / n
                funcArr[i] = func(y * temp1 + temp2);
            }
        }
        const T factor(2.0 * n_1);
        for (size_t i = 0; i < n; ++i) {
            T sum(0);
            for (size_t j = 0; j < n; ++j)
                sum += funcArr[j] * std::cos(pi * i * (j + 0.5) * n_1);
            coeff[i] = factor * sum;
        }
    }
    /**
     * Fit a even function, the performance is better than chebyshev_fit()
     * 
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:143
     */
    template<Scalar T, Vector U, class Function>
    void chebyshev_fit_even(const T& from, const T& to, LValueVector<U>& coeff, Function func) {
        constexpr double pi = M_PI;
        assert(from < to);

        const size_t n = coeff.getLength();
        const size_t n_2 = n << 1U;
        const double n_2_1 = 1.0 / n_2;
        VectorND<T> funcArr(n_2);
        /* Fill func arr */ {
            const T& temp1 = (to - from) * T(0.5);
            const T& temp2 = (to + from) * T(0.5);
            for (size_t i = 0; i < n_2; ++i) {
                const T y = T(std::cos(pi * (i + 0.5) * n_2_1));
                funcArr[i] = func(y * temp1 + temp2);
            }
        }
        const T factor(4.0 * n_2_1);
        for (size_t i = 0; i < n; ++i) {
            T sum(0);
            for (size_t j = 0; j < n; ++j)
                sum += funcArr[j] * std::cos(2 * pi * i * (j + 0.5) * n_2_1);
            coeff[i] = factor * sum;
        }
    }
    /**
     * Fit a odd function, the performance is better than chebyshev_fit()
     * 
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:143
     */
    template<Scalar T, Vector U, class Function>
    void chebyshev_fit_odd(const T& from, const T& to, LValueVector<U>& coeff, Function func) {
        constexpr double pi = M_PI;
        assert(from < to);

        const size_t n = coeff.getLength();
        const size_t n_2 = n << 1U;
        const double n_2_1 = 1.0 / n_2;
        VectorND<T> funcArr(n_2);
        /* Fill func arr */ {
            const T& temp1 = (to - from) * T(0.5);
            const T& temp2 = (to + from) * T(0.5);
            for (size_t i = 0; i < n_2; ++i) {
                const T y = T(std::cos(pi * (i + 0.5) * n_2_1));
                const T x = y * temp1 + temp2;
                funcArr[i] = func(x) / x;
            }
        }
        const T factor(4.0 * n_2_1);
        for (size_t i = 0; i < n; ++i) {
            T sum(0);
            for (size_t j = 0; j < n; ++j)
                sum += funcArr[j] * std::cos(2 * pi * i * (j + 0.5) * n_2_1);
            coeff[i] = factor * sum;
        }
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:142
     */
    template<Scalar T, Vector U>
    T chebyshev_calc(const T& from, const T& to, const U& coeff, const T& x) {
        assert(from <= x && x <= to);
        const T y = T(x * T(2) - to - from) / T(to - from);
        const T y2 = y * T(2);
        T d1(0), d2(0);
        for (size_t i = coeff.getLength() - 1; i > 0; --i) {
            const T temp = d1;
            d1 = y2 * d1 - d2 + coeff.calc(i);
            d2 = temp;
        }
        return y * d1 - d2 + coeff.calc(0) * T(0.5);
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:143
     */
    template<Scalar T, Vector U>
    T chebyshev_calc_even(const T& from, const T& to, const U& coeff, const T& x) {
        assert(from <= x && x <= to);
        return chebyshev_calc<T, U>(from, to, coeff, square(x) * T(2) - T(1));
    }

    template<Scalar T, Vector U>
    T chebyshev_calc_odd(const T& from, const T& to, const U& coeff, const T& x) {
        assert(from <= x && x <= to);
        return chebyshev_calc_even<T, U>(from, to, coeff, x) * x;
    }
}
