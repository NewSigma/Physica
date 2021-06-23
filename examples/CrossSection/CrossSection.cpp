/*
 * Copyright 2020 WeiBo He.
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
#include <QtWidgets/QApplication>
#include <Physica/Core/Math/Calculus/ODE/ODESolver.h>
#include <Physica/Core/Physics/PhyConst.h>
#include <Physica/Gui/Plot/Plot.h>

using namespace Physica::Core;
using namespace Physica::Core::Physics;
using namespace Physica::Gui;
using T = Scalar<Double, false>;
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013.3
 */
int main(int argc, char** argv) {
    using ODE = ODESolver<T, Vector<T>>;
    const T m2_h_2 = 2 * PhyConst::electroMass / (PhyConst::reducedPlanck * PhyConst::reducedPlanck);
    T stepSize(0.01);
    T rmax = 1;
    T energy = 1;
    const T wave_length = sqrt(T(2) / (T(PhyConst::electroMass) * energy)) * T(M_PI * PhyConst::reducedPlanck);
    constexpr double radialNum = 0;

    ODE solver(0.001, rmax + wave_length, stepSize, {0});
    solver.degenerate_numerov([&](T r) -> T {
        const T V = r > rmax ? 0 : -5;
        return m2_h_2 * (V - energy) + T(radialNum * (radialNum + 1.0)) / square(r);
    }, {0.01});
    const auto& x = solver.getX();
    const auto& solution = solver.getSolution();

    Vector<T> waveFunc(solution.getColumn());
    for (size_t i = 0; i < solution.getColumn(); ++i)
        waveFunc[i] = solution[i][0];

    QApplication app(argc, argv);
    Plot* plot = new Plot(x, waveFunc);
    plot->show();
    return QApplication::exec();
}
