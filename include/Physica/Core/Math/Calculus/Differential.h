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

namespace Physica::Core {
    //Reference: Numerical Recipes in C++
    template<class ScalarType>
    class Differential {
    public:
        template<class Function>
        [[nodiscard]] static ScalarType forward(Function func, const ScalarType& x, const ScalarType& step);
        template<class Function>
        [[nodiscard]] static ScalarType backward(Function func, const ScalarType& x, const ScalarType& step);
        template<class Function>
        [[nodiscard]] static ScalarType doublePoint(Function func, const ScalarType& x, const ScalarType& step);
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
}
