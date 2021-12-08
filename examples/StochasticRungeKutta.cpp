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
#include <random>
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

int main(int argc, char** argv) {
    using T = Scalar<Double, false>;
    using ODE = ODESolver<T, Vector<T>>;
    constexpr double stepSize = 0.01;
    constexpr double D = 1;
    constexpr double gamma = 1;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution phi{};

    ODE solver(0, 50, stepSize, {0, 0});
    solver.stochasticRungeKutta2([](T x, const Vector<T>& y) -> Vector<T> { (void)x; return {y[1], -T(gamma) * y[1]}; },
                                [&](T x, const Vector<T>& y) -> Vector<T> { (void)x; (void)y; return {0, sqrt(2 * stepSize * D) * phi(gen)}; });
    const auto& t = solver.getX();
    const auto& solution = solver.getSolution();

    const size_t length = solution.getColumn();
    Vector<T> x, v;
    x.resize(length);
    v.resize(length);
    for (size_t i = 0; i < length; ++i) {
        x[i] = solution(0, i);
        v[i] = solution(1, i);
    }

    QApplication app(argc, argv);
    /* t-x */ {
        Plot* plot = new Plot();
        plot->spline(t, x);
        auto& chart = *plot->chart();
        chart.legend()->hide();
        chart.setTitle("t_x");
        chart.createDefaultAxes();
        auto axes = chart.axes();
        axes[0]->setTitleText("t");
        axes[1]->setTitleText("x");
        plot->show();    
    }
    /* t-v */ {
        Plot* plot = new Plot();
        plot->spline(t, v);
        auto& chart = *plot->chart();
        chart.legend()->hide();
        chart.setTitle("t_v");
        chart.createDefaultAxes();
        auto axes = chart.axes();
        axes[0]->setTitleText("t");
        axes[1]->setTitleText("v");
        plot->show();
    }
    /* x-v */ {
        Plot* plot = new Plot();
        plot->scatter(x, v);
        auto& chart = *plot->chart();
        auto* series = static_cast<QScatterSeries*>(chart.series().first());
        series->setMarkerSize(5);
        chart.legend()->hide();
        chart.setTitle("x_v");
        chart.createDefaultAxes();
        auto axes = chart.axes();
        axes[0]->setTitleText("x");
        axes[1]->setTitleText("v");
        plot->show();   
    }
    return QApplication::exec();
}