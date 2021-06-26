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
#include "Physica/Core/Math/Algebra/BoolAlgebra/BoolMatrix.h"
#include "Physica/Core/Math/Statistics/LinearFit.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

using T = Scalar<Double, false>;
using ODE = ODESolver<T, Vector<T>>;

constexpr double side = 4;

size_t fillPhase(const Vector<T>& x, const Vector<T>& p, int dim) {
    constexpr double x0 = -2;
    constexpr double y0 = -2;

    double deltaWidth = side / dim;
    double deltaHeight = side / dim;
    size_t fillCount = 0;
    BoolMatrix boolMat(dim, dim, false);

    for (size_t i = 0; i < x.getLength(); ++i) {
        int m = floor((x[i].getTrivial() - x0) / deltaWidth);
        int n = floor((p[i].getTrivial() - y0) / deltaHeight);
        if (0 <= m && m < dim) {
            if (0 <= n && n < dim) {
                fillCount += !boolMat[n][m];
                boolMat.setValue(n, m, true);
            }
        }
    }
    return fillCount;
}
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:11-12
 */
int main(int argc, char** argv) {
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

    /* Get fractal dimention */ {
        constexpr unsigned int i_start = 2U;
        constexpr unsigned int i_end = 10U;
        Vector<T> x1(i_end - i_start);
        Vector<T> y1(i_end - i_start);
        for (unsigned int i = i_start; i < i_end; ++i) {
            x1[i - i_start] = log(side / (1U << i));
            y1[i - i_start] = log(fillPhase(x, p, 1U << i));
        }
        auto fit = linearFit(x1, y1);
        std::cout << "Fractal dimention: " << -fit.first.getTrivial() << std::endl;
    }


    QApplication app(argc, argv);
    Plot* t_x = new Plot();
    t_x->spline(t, x);
    t_x->show();
    Plot* t_p = new Plot();
    t_p->spline(t, p);
    t_p->show();
    Plot* x_p = new Plot();
    x_p->spline(x, p);
    x_p->show();
    return QApplication::exec();
}