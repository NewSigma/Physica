/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Calculus/ODE/SRK2.h"
#include "Physica/Core/Math/Statistics/LinearFit.h"

using namespace Physica;
using RandomSource = Random<MT19937, 16507335086333683534UL>;
/**
 * We solve Eq. 5.1 of [1] and the slope of the result is expected to be lambda.
 *
 * Reference:
 * [1] Phys. Rev. A 45, 600 (1992); https://doi.org/10.1103/PhysRevA.45.600
 */
int main() {
    using T = float64;
    using XVector = Vector1D<T>;
    constexpr double stepSize = 0.1;
    constexpr double diffuseD = 0.1;
    constexpr double lambda = 1;
    constexpr double t_max = 4;
    constexpr size_t iteration = 10000;

    SRK2<T, 1> solver(0, t_max, stepSize, {1});
    VectorND<T> y(solver.getNumStep(), 0);
    for (size_t i = 0; i < iteration; ++i) {
        solver.solve([](T, const XVector& x) -> XVector { return {-T(lambda) * x[0]}; },
                     [](T, const XVector&) -> XVector {
                         return {T(lambda) * sqrt(2 * stepSize * diffuseD) * T::random_normal<RandomSource>()};
                     });
        y.toNextMean(i, solver.getSolution().row(0));
    }

    const VectorND<T> t = VectorND<T>::linspace(0, t_max, y.getLength());
    const VectorND<T> ln_y = ln(y);
    const auto param = LinearFit<T>::fit(t, ln_y);
    const auto devia = LinearFit<T>::deviation(t, ln_y, param);
    if (abs(-param.first - lambda) > T(3) * devia.first)
        return 1;
    return 0;
}
