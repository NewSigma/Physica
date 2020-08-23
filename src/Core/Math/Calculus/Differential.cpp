/*
 * Copyright 2020 WeiBo He.
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
#include "Physica/Core/Math/Calculus/Differential.h"

namespace Physica::Core {
    Differential::Differential(TreeFunction<> func, MultiScalar at, MultiScalar stepSize)
            : func(std::move(func)), at(std::move(at)), stepSize(std::move(stepSize)) {}
    /*!
     * Optimize: if \at is much larger than \stepsize, the result will be 0. May be use talor series
     * and expend the function to the first order.
     */
    MultiScalar Differential::solve(DifferentialMethod method) {
        MultiScalar result;
        switch(method) {
            case DoublePoint:
                result = (func(at + stepSize) - func(at - stepSize)) / (stepSize << 1);
                break;
            case Forward:
                result = (func(at + stepSize) - func(at)) / stepSize;
                break;
            case Backward:
                result = (func(at) - func(at - stepSize)) / stepSize;
                break;
        }
        return result;
    }
}