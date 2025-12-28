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
#include <iostream>
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"
#include "Physica/Core/Physics/Experiment/DimEstimator.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using T = float64;
using ODE = ODESolver<T, 2>;
/**
 * Reference:
 * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:11-12
 */
int main(int argc, char** argv) {
    ODE solver(0, 200, 0.01, {0.5, 0});
    solver.rungeKutta4([](const T& t, const Vector2D<T>& x) { return Vector2D<T>{x[1], -x[1] * T(0.1) + (T(0.5) - T(2) * x[0] * x[0]) * x[0] + T(2) * cos(T(2.4) * t)}; });
    const auto& t = solver.getX();
    const auto& solution = solver.getSolution();

    VectorND<T> x{};
    x.resize(solution.getCol());
    VectorND<T> p{};
    p.resize(solution.getCol());

    for (size_t i = 0; i < solution.getCol(); ++i) {
        x[i] = solution[0, i];
        p[i] = solution[1, i];
    }

    /* Get fractal dimention */ {
        const size_t length = 32;
        using VectorType = VectorND<T>;
        DenseMatrix<T> trans = solution.transpose();
        const VectorType r = exp(VectorType::linspace(-4, -1, length));
        const T dim = DimEstimator::corrDimen(trans, r);
        std::cout << "Effective dimention: " << dim << '\n';
    }

    QApplication app(argc, argv);
    /* t-x */ {
        Plot* t_x = new Plot(0, 200, -1.5, 1.5, 50, 1);
        t_x->spline(t, x);
        t_x->getChart()->legend()->hide();
        auto* axisX = t_x->getAxisX();
        auto* axisY = t_x->getAxisY();
        axisX->setTitleText("t");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("x");
        axisY->setLabelFormat("%d");
        t_x->show();
    }
    /* t-p */ {
        Plot* t_p = new Plot(0, 200, -2.5, 2, 50, 1);
        t_p->spline(t, p);
        t_p->getChart()->legend()->hide();
        auto* axisX = t_p->getAxisX();
        auto* axisY = t_p->getAxisY();
        axisX->setTitleText("t");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("p");
        axisY->setLabelFormat("%d");
        t_p->show();
    }
    /* x-p */ {
        Plot* x_p = new Plot(-1.5, 1.5, -2.5, 2, 1, 1);
        x_p->spline(x, p);
        x_p->getChart()->legend()->hide();
        auto* axisX = x_p->getAxisX();
        auto* axisY = x_p->getAxisY();
        axisX->setTitleText("x");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("p");
        axisY->setLabelFormat("%d");
        x_p->show();
    }
    return QApplication::exec();
}