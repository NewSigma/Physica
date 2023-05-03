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
#include <iostream>
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"

using namespace Physica::Core;

int main() {
    using T = Scalar<Double, false>;
    constexpr double stepSize = 0.0001;
    /**
     * y' = y
     * y[0] = 1
     */
    {
        ODESolver<T, 1> solver(0, 3, stepSize, {1});
        solver.rungeKutta4([](T x, const Vector<T, 1>& y) -> Vector<T, 1> { (void)x; return y; });
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getColumn(); ++i) {
            const auto& solVector = solution[i];
            const auto answer = exp(x[i].getTrivial());
            if (fabs(solVector[0].getTrivial() - answer) > 1E-11 * fabs(answer))
                return 1;
        }
    }
    /**
     * y'' + y = 0
     * y[0] = 0  y'[0] = 1
     */
    {
        ODESolver<T, 2> solver(0, 3, stepSize, {0, 1});
        solver.rungeKutta4([](T x, const Vector<T, 2>& y) -> Vector<T, 2> { (void)x; return {y[1], -y[0]}; });
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getColumn(); ++i) {
            const auto& x_i = x[i];
            const auto& solVector = solution[i];
            const auto answer1 = sin(x_i.getTrivial());
            if (fabs(solVector[0].getTrivial() - answer1) > 1E-10 * fabs(answer1))
                return 1;
            const auto answer2 = cos(x_i.getTrivial());
            if (fabs(solVector[1].getTrivial() - answer2) > 1E-7 * fabs(answer2))
                return 1;
        }
    }
    /**
     * y'' = y
     * y[0] = 1
     */
    {
        ODESolver<T, 2> solver(0, 3, stepSize, {1});
        solver.verlet([](T x, const T& y) -> T { (void)x; return y; }, {exp(stepSize)});
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getColumn(); ++i) {
            const auto& solVector = solution[i];
            const auto answer = exp(x[i].getTrivial());
            if (fabs(solVector[0].getTrivial() - answer) > 1E-9 * fabs(answer))
                return 1;
        }
    }
    /**
     * y'' + y = 0
     * y[0] = 0  y'[0] = 1
     */
    {
        ODESolver<T, 2> solver(0, 3, stepSize, {0});
        solver.degenerate_numerov([](T x) -> T { (void)x; return -1; }, {1});
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getColumn(); ++i) {
            const auto& x_i = x[i];
            const auto& solVector = solution[i];
            const auto answer = sin(x_i.getTrivial());
            if (fabs(solVector[0].getTrivial() - answer) > 1E-9 * fabs(answer))
                return 1;
        }
    }
    /**
     * y'' = (x^2 - 3)y
     * y[0] = 0  y'[0] = 1
     * Comes from Check 1 in [1]
     * 
     * Reference:
     * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:20
     */
    {
        ODESolver<T, 2> solver(0, 2, stepSize, {0});
        solver.degenerate_numerov([](T x) -> T { return square(x) - T(3); }, {1});
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getColumn(); ++i) {
            const auto& x_i = x[i];
            const auto& solVector = solution[i];
            const auto answer = (x_i * exp(-square(x_i) / 2)).getTrivial();
            if (fabs(solVector[0].getTrivial() - answer) > 1E-8 * fabs(answer))
                return 1;
        }
    }
    return 0;
}