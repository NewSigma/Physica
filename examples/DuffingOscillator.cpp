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
#include "Physica/Core/Physics/Experiment/DimEstimator.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

using T = Scalar<Double, false>;
using ODE = ODESolver<T>;
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
        const size_t length = 32;
        using VectorType = Vector<T>;
        DenseMatrix<T> trans = solution.transpose();
        const VectorType r = exp(VectorType::linspace(-4, -1, length));
        const T dim = DimEstimator::corrDimen(trans, r);
        std::cout << "Effective dimention: " << dim << std::endl;
    }

    QApplication app(argc, argv);
    /* t-x */ {
        Plot* t_x = new Plot();
        t_x->spline(t, x);
        auto& chart = *t_x->chart();
        chart.setTitle("t-x");
        chart.legend()->hide();
        chart.createDefaultAxes();
        chart.axes(Qt::Horizontal).first()->setTitleText("t");
        chart.axes(Qt::Vertical).first()->setTitleText("x");
        t_x->show();
    }
    /* t-p */ {
        Plot* t_p = new Plot();
        t_p->spline(t, p);
        auto& chart = *t_p->chart();
        chart.setTitle("t-p");
        chart.legend()->hide();
        chart.createDefaultAxes();
        chart.axes(Qt::Horizontal).first()->setTitleText("t");
        chart.axes(Qt::Vertical).first()->setTitleText("p");
        t_p->show();
    }
    /* x-p */ {
        Plot* x_p = new Plot();
        x_p->spline(x, p);
        auto& chart = *x_p->chart();
        chart.setTitle("x-p");
        chart.legend()->hide();
        chart.createDefaultAxes();
        chart.axes(Qt::Horizontal).first()->setTitleText("x");
        chart.axes(Qt::Vertical).first()->setTitleText("p");
        x_p->show();
    }
    return QApplication::exec();
}