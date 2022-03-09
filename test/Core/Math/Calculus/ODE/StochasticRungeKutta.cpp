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
#include <random>
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"

using namespace Physica::Core;
/**
 * Reference:
 * [1] R. L. Honeycutt, Stochastic Runge-Kutta algorithm: I. White noise, Phys. Rev. A 45, 600 (1992).
 */
int main() {
    using ScalarType = Scalar<Double, false>;
    using ODE = ODESolver<ScalarType>;
    using XVector = Vector<ScalarType, 1>;
    constexpr double stepSize = 0.1;
    constexpr double D = 0.1;
    constexpr double lambda = 1;
    constexpr double t_max = 4;
    constexpr size_t iteration = 5000;

    const size_t count = t_max / stepSize;
    std::random_device rd{};
    Vector<ScalarType> x(count, 0);
    for (size_t i = 0; i < iteration; ++i) {
        std::mt19937 gen{i};
        std::normal_distribution phi{};

        ODE solver(0, t_max, stepSize, {1});
        solver.stochasticRungeKutta2([](ScalarType t, const XVector& x) -> XVector { (void)t; return {-ScalarType(lambda) * x[0]}; },
                                    [&](ScalarType t, const XVector& x) -> XVector { (void)t; (void)x; return {ScalarType(lambda) * sqrt(2 * stepSize * D) * phi(gen)}; });
        x += solver.getSolution().row(0);
    }
    x *= reciprocal(ScalarType(iteration));
    const Vector<ScalarType> log_x = ln(x);
    const Vector<ScalarType> t = Vector<ScalarType>::linspace(0, t_max - stepSize, count);
    const Vector<ScalarType> log_x_theory = ScalarType(-lambda) * t;
    if (abs(log_x - log_x_theory).max() > ScalarType(0.178))
        return 1;
    return 0;
}
