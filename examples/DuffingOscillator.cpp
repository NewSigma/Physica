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
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013.3
 */
int main(int argc, char** argv) {
    using namespace Physica::Core;
    using namespace Physica::Gui;

    using T = Scalar<Double, false>;
    using ODE = ODESolver<T, Vector<T>>;

    ODE solver(0, 200, 0.01, {0.5, 0});
    solver.rungeKutta4([](const T& t, const Vector<T>& x) { return Vector<T>{x[1], -x[1] * T(0.1) + (T(0.5) - T(2) * x[0] * x[0]) * x[0] + T(2) * cos(T(2.4) * t)}; });
    const auto& t = solver.getX();
    const auto& solution = solver.getSolution();
    Vector<T> x{};
    x.resize(solution.getColumn());
    Vector<T> p{};
    p.resize(solution.getColumn());
    for (size_t i = 0; i < solution.getColumn(); ++i) {
        x[i] = solution(0, i);
        p[i] = solution(1, i);
    }

    QApplication app(argc, argv);
    Plot* t_x = new Plot(t, x);
    t_x->show();
    Plot* t_p = new Plot(t, p);
    t_p->show();
    Plot* x_p = new Plot(x, p);
    x_p->show();
    return QApplication::exec();
}