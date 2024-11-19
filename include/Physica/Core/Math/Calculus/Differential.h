/*
 * Copyright 2020-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h>

namespace Physica::Core {
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:137-140
     */
    template<Scalar T>
    class Differential {
    public:
        template<class Function>
        [[nodiscard]] static T forward(Function func, T x, T step);
        template<class Function>
        [[nodiscard]] static T backward(Function func, T x, T step);
        template<class Function>
        [[nodiscard]] static T doublePoint(Function func, T x, T step);
        template<class Function>
        [[nodiscard]] static T ridders(Function func, T x, T step);
    };

    template<Scalar T>
    template<class Function>
    T Differential<T>::forward(Function func, T x, T step) {
        return (func(x + step) - func(x)) / step;
    }

    template<Scalar T>
    template<class Function>
    T Differential<T>::backward(Function func, T x, T step) {
        return (func(x) - func(x - step)) / step;
    }

    template<Scalar T>
    template<class Function>
    T Differential<T>::doublePoint(Function func, T x, T step) {
        return (func(x + step) - func(x - step)) / (step  * T(2));
    }

    template<Scalar T>
    template<class Function>
    T Differential<T>::ridders(Function func, T x, T step) {
        constexpr static size_t TableSize = 10;
        constexpr static double Factor = 1.4;
        constexpr static double RepFactor = 1 / Factor;
        constexpr static double Factor2 = Factor * Factor;
        constexpr static double Tolerance = 2;

        assert(step.isPositive());
        DenseSymmMatrix<T, TableSize> table(TableSize);
        table(0, 0) = doublePoint(func, x, step);
        T step_now = step;
        T error = std::numeric_limits<T>::max();
        T result{};
        
        for (size_t i = 1; i < TableSize; ++i) {
            step_now *= T(RepFactor);
            table(0, i) = doublePoint(func, x, step_now);

            T factor2 = Factor2;
            for (size_t j = 1; j <= i; ++j) {
                table(j, i) = (table(j - 1, i) * factor2 - table(j - 1, i - 1)) / (factor2 - 1);
                factor2 *= T(Factor2);
                const T error_now = std::max(abs(table(j, i) - table(j - 1, i)), abs(table(j, i) - table(j - 1, i - 1)));
                if (error_now <= error) {
                    error = error_now;
                    result = table(j, i);
                }
            }
            if (abs(table(i, i) - table(i - 1, i - 1)) >= error * Tolerance)
                break;
        }
        return result;
    }
}
