/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Math/Calculus/Integrate/Vegas.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using namespace Physica::Gui;
using RandomType = Random<MT19937>;
using T = float64;

T func(const VectorND<T>& x) {
    return exp(x.squaredNorm() * T(-0.5));
}

void plotCompressRate() {
    constexpr double r = 3;
    const VectorND<T> from(8, T(-r));
    const VectorND<T> to(8, T(r));

    Plot* plot = new Plot(0, 100, -5, -0.5, 25, 1);
    plot->getLegend().setAlignment(Qt::AlignRight);
    plot->getLegend().show();
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%d");
    plot->getAxisX()->setTitleText("NumRefine");
    plot->getAxisY()->setTitleText("ln(L)");

    double rates[3]{0.1, 0.2, 0.5};
    const char* names[3]{"0.1", "0.2", "0.5"};
    for (int i = 0; i < 3; ++i) {
        Vegas<T, false> vegas(from, to, 100, 100000, 10, rates[i]);
        vegas.integral<decltype(func), RandomType, ThreadExecutor>(func);
        VectorND<T> vars = ln(vegas.getLoss());
        plot->line(vars).setName(names[i]);
    }
    plot->show();
    QApplication::exec();
}

void plotNumPoint() {
    constexpr double r = 3;
    const VectorND<T> from(8, T(-r));
    const VectorND<T> to(8, T(r));

    Plot* plot = new Plot(0, 1000, -9, -1, 250, 2);
    plot->getLegend().setAlignment(Qt::AlignRight);
    plot->getLegend().show();
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%d");
    plot->getAxisX()->setTitleText("NumRefine");
    plot->getAxisY()->setTitleText("ln(L)");

    int points[3]{10, 100, 1000};
    const char* names[3]{"10", "100", "1000"};
    for (int i = 0; i < 3; ++i) {
        Vegas<T, false> vegas(from, to, 1000, 100000, points[i], 0.1);
        vegas.integral<decltype(func), RandomType, ThreadExecutor>(func);
        VectorND<T> vars = ln(vegas.getLoss());
        plot->line(vars).setName(names[i]);
    }
    plot->show();
    QApplication::exec();
}

void plotNumSample() {
    constexpr double r = 3;
    const VectorND<T> from(8, T(-r));
    const VectorND<T> to(8, T(r));

    Plot* plot = new Plot(0, 1000, -9, -5, 250, 1);
    plot->getLegend().setAlignment(Qt::AlignRight);
    plot->getLegend().show();
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%d");
    plot->getAxisX()->setTitleText("NumRefine");
    plot->getAxisY()->setTitleText("ln(L)");

    int samples[2]{10000, 100000};
    const char* names[2]{"10<sup>4</sup>", "10<sup>5</sup>"};
    for (int i = 0; i < 2; ++i) {
        Vegas<T, false> vegas(from, to, 1000, samples[i], 1000, 0.1);
        vegas.integral<decltype(func), RandomType, ThreadExecutor>(func);
        VectorND<T> vars = ln(vegas.getLoss());
        plot->line(vars).setName(names[i]);
    }
    plot->show();
    QApplication::exec();
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    QApplication app(argc, argv);
    plotCompressRate();
    plotNumPoint();
    plotNumSample();
    return 0;
}
