/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"
#include "Test.h"

using namespace Physica;

int main() {
    using T = float64;
    constexpr double stepSize = 0.0001;
    /**
     * y' = y
     * y[0] = 1
     */
    {
        ODESolver<T, 1> solver(0, 3, stepSize, {1});
        solver.rungeKutta4([](T x, const Vector1D<T>& y) -> Vector1D<T> { (void)x; return y; });
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getCol(); ++i) {
            const T result = solution.col(i)[0];
            const T answer = exp(x[i]);
            expect(scalarNear(result, answer, 1E-11));
        }
    }
    /**
     * y'' + y = 0
     * y[0] = 0  y'[0] = 1
     */
    {
        ODESolver<T, 2> solver(0, 3, stepSize, {0, 1});
        solver.rungeKutta4([](T x, const Vector2D<T>& y) -> Vector2D<T> { (void)x; return {y[1], -y[0]}; });
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getCol(); ++i) {
            const auto& x_i = x[i];
            const auto solVector = solution.col(i);
            const auto answer1 = sin(x_i);
            expect(scalarNear(solVector[0], answer1, 1E-10));

            const auto answer2 = cos(x_i);
            expect(scalarNear(solVector[1], answer2, 1E-7));
        }
    }
    /**
     * y'' = y
     * y[0] = 1
     */
    {
        ODESolver<T, 1> solver(0, 3, stepSize, {1});
        solver.verlet([](T x, const T& y) -> T { (void)x; return y; }, {exp(stepSize)});
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getCol(); ++i) {
            const auto solVector = solution.col(i);
            const auto answer = exp(x[i]);
            expect(scalarNear(solVector[0], answer, 1E-8));
        }
    }
    /**
     * y'' + y = 0
     * y[0] = 0  y'[0] = 1
     */
    {
        ODESolver<T, 2> solver(0, 3, stepSize, {0, 1});
        solver.degenerate_numerov([](T x) -> T { (void)x; return -1; }, {1});
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getCol(); ++i) {
            const auto& x_i = x[i];
            const auto solVector = solution.col(i);
            const auto answer = sin(x_i);
            expect(scalarNear(solVector[0], answer, 1E-8));
        }
    }
    /**
     * y'' = (x^2 - 3)y
     * y[0] = 0  y'[0] = 1
     * Comes from Check 1 in [1]
     *
     * Reference:
     * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:20
     */
    {
        ODESolver<T, 2> solver(0, 2, stepSize, {0, 1});
        solver.degenerate_numerov([](T x) -> T { return square(x) - T(3); }, {1});
        const auto& x = solver.getX();
        const auto& solution = solver.getSolution();
        for (size_t i = 0; i < solution.getCol(); ++i) {
            const auto& x_i = x[i];
            const auto solVector = solution.col(i);
            const auto answer = x_i * exp(-square(x_i) / 2);
            expect(scalarNear(solVector[0], answer, 1E-8));
        }
    }
    return 0;
}
