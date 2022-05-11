/*
 * Copyright 2020-2021 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009:137-140
     */
    template<class ScalarType>
    class Differential {
    public:
        template<class Function>
        [[nodiscard]] static ScalarType forward(Function func, const ScalarType& x, const ScalarType& step);
        template<class Function>
        [[nodiscard]] static ScalarType backward(Function func, const ScalarType& x, const ScalarType& step);
        template<class Function>
        [[nodiscard]] static ScalarType doublePoint(Function func, const ScalarType& x, const ScalarType& step);
        template<class Function>
        [[nodiscard]] static ScalarType ridders(Function func, const ScalarType& x, const ScalarType& step);
    };

    template<class ScalarType>
    template<class Function>
    ScalarType Differential<ScalarType>::forward(Function func, const ScalarType& x, const ScalarType& step) {
        return (func(x + step) - func(x)) / step;
    }

    template<class ScalarType>
    template<class Function>
    ScalarType Differential<ScalarType>::backward(Function func, const ScalarType& x, const ScalarType& step) {
        return (func(x) - func(x - step)) / step;
    }

    template<class ScalarType>
    template<class Function>
    ScalarType Differential<ScalarType>::doublePoint(Function func, const ScalarType& x, const ScalarType& step) {
        return (func(x + step) - func(x - step)) / (step  * ScalarType::Two());
    }

    template<class ScalarType>
    template<class Function>
    ScalarType Differential<ScalarType>::ridders(Function func, const ScalarType& x, const ScalarType& step) {
        constexpr static size_t TableSize = 10;
        constexpr static double Factor = 1.4;
        constexpr static double RepFactor = 1 / Factor;
        constexpr static double Factor2 = Factor * Factor;
        constexpr static double Tolerance = 2;

        assert(step.isPositive());
        DenseSymmMatrix<ScalarType, TableSize, TableSize> table(TableSize);
        table(0, 0) = doublePoint(func, x, step);
        ScalarType step_now = step;
        ScalarType error = std::numeric_limits<ScalarType>::max();
        ScalarType result{};
        
        for (size_t i = 1; i < TableSize; ++i) {
            step_now *= ScalarType(RepFactor);
            table(0, i) = doublePoint(func, x, step_now);

            ScalarType factor2 = Factor2;
            for (size_t j = 1; j <= i; ++j) {
                table(j, i) = (table(j - 1, i) * factor2 - table(j - 1, i - 1)) / (factor2 - 1);
                factor2 *= ScalarType(Factor2);
                const ScalarType error_now = std::max(abs(table(j, i) - table(j - 1, i)), abs(table(j, i) - table(j - 1, i - 1)));
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
