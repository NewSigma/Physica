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
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

using T = Scalar<Double, false>;
using ODE = ODESolver<T, Vector<T>>;

int main(int argc, char** argv) {
    const T radius = 0.005;
    const T rmin = 0.000001;
    const T stepSize = rmin;
    const T sigma = 0.074;
    const T rho = 1000;
    const T g = 9.8;
    const T lambda = -0.0300034492;
    ODE solver(-radius, -rmin, stepSize, {0, 1});
    solver.rungeKutta4([&](const T& r, const Vector<T>& v) {
            const T momentum = v[1];
            const T momentum_2_1 = momentum * momentum + T(1);
            const T sqrt_momentum_2_1 = sqrt(momentum_2_1);
            const T term1 = (momentum_2_1 * sqrt_momentum_2_1 * (rho * g * v[0] + lambda * rho)) / sigma;
            const T term2 = momentum_2_1 * momentum / r;
            return Vector<T>{momentum, term1 - term2};
        });
    const auto& r = solver.getX();
    const auto& solution = solver.getSolution();

    const size_t length = solution.getColumn();
    Vector<T> z{};
    z.resize(length);

    T volumeHelper(0);
    T tangent_min = abs(solution(1, 0));
    for (size_t i = 0; i < length; ++i) {
        T temp1 = solution(0, i);
        z[i] = temp1;
        volumeHelper += temp1 * r[i];
        T temp2 = abs(solution(1, i));
        tangent_min = temp2 < tangent_min ? temp2 : tangent_min;
    }
    std::cout << "Volume is " << (-volumeHelper * stepSize) << " m^3" << std::endl;
    std::cout << "Minimum tangent is " << tangent_min << std::endl;

    QApplication app(argc, argv);
    Plot* r_z = new Plot(r, z);
    r_z->show();
    return QApplication::exec();
}
