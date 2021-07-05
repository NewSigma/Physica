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

int main(int argc, char** argv) {
    using T = Scalar<Double, false>;
    using ODE = ODESolver<T, Vector<T>>;
    constexpr double stepSize = 0.0001;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    ODE solver(0, 2, stepSize, {0});
    solver.stochasticRungeKutta2([](T x, const Vector<T>& y) -> Vector<T> { (void)x; return -y; }, gen, 1);
    const auto& x = solver.getX();
    const auto& solution = solver.getSolution();

    Vector<T> v{};
    v.resize(solution.getColumn());
    for (size_t i = 0; i < solution.getColumn(); ++i) {
        v[i] = solution(0, i);
    }

    QApplication app(argc, argv);
    Plot* t_x = new Plot();
    t_x->spline(x, v);
    t_x->show();
    return QApplication::exec();
}