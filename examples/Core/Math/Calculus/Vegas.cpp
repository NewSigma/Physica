/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Calculus/Integrate/Vegas.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Algorithm/Thread.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using RandomSource = Random<>;
using T = float64;
using Cube = Vegas<T, false>::Cube;

namespace {
    T func(const VectorND<T>& x) {
        return exp(x.squaredNorm() * T(-0.5));
    }

    void plotCompressRate() {
        constexpr double r = 3;
        const Cube cube{VectorND<T>(8, T(-r)), VectorND<T>(8, T(r))};

        Plot* plot = new Plot(0, 100, -5, -0.5, 25, 1);
        plot->getLegend().setAlignment(Qt::AlignRight);
        plot->getLegend().show();
        plot->getAxisX()->setLabelFormat("%d");
        plot->getAxisY()->setLabelFormat("%d");
        plot->getAxisX()->setTitleText("NumRefine");
        plot->getAxisY()->setTitleText("ln(L)");

        Array<double, 3> rates{0.1, 0.2, 0.5};
        Array<const char*, 3> names{"0.1", "0.2", "0.5"};
        for (int i = 0; i < 3; ++i) {
            Vegas<T, false> vegas(cube, 100000, 10, rates[i]);
            const auto process = vegas.integral<RandomSource, Thread>(func, 100);
            VectorND<T> vars = ln(process.getLoss());
            plot->line(vars).setName(names[i]);
        }
        plot->show();
        QApplication::exec();
    }

    void plotNumPoint() {
        constexpr double r = 3;
        const Cube cube{VectorND<T>(8, T(-r)), VectorND<T>(8, T(r))};

        Plot* plot = new Plot(0, 1000, -9, -1, 250, 2);
        plot->getLegend().setAlignment(Qt::AlignRight);
        plot->getLegend().show();
        plot->getAxisX()->setLabelFormat("%d");
        plot->getAxisY()->setLabelFormat("%d");
        plot->getAxisX()->setTitleText("NumRefine");
        plot->getAxisY()->setTitleText("ln(L)");

        Array<int, 3> points{10, 100, 1000};
        for (int point : points) {
            Vegas<T, false> vegas(cube, 100000, point, 0.1);
            const auto process = vegas.integral<RandomSource, Thread>(func, 1000);
            VectorND<T> vars = ln(process.getLoss());
            plot->line(vars).setName(std::format("{}", point).c_str());
        }
        plot->show();
        QApplication::exec();
    }

    void plotNumSample() {
        constexpr double r = 3;
        const Cube cube{VectorND<T>(8, T(-r)), VectorND<T>(8, T(r))};

        Plot* plot = new Plot(0, 1000, -9, -5, 250, 1);
        plot->getLegend().setAlignment(Qt::AlignRight);
        plot->getLegend().show();
        plot->getAxisX()->setLabelFormat("%d");
        plot->getAxisY()->setLabelFormat("%d");
        plot->getAxisX()->setTitleText("NumRefine");
        plot->getAxisY()->setTitleText("ln(L)");

        Array<int, 2> samples{10000, 100000};
        Array<const char*, 2> names{"10<sup>4</sup>", "10<sup>5</sup>"};
        for (int i = 0; i < 2; ++i) {
            Vegas<T, false> vegas(cube, samples[i], 1000, 0.1);
            const auto process = vegas.integral<RandomSource, Thread>(func, 1000);
            VectorND<T> vars = ln(process.getLoss());
            plot->line(vars).setName(names[i]);
        }
        plot->show();
        QApplication::exec();
    }
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    QApplication app(argc, argv);
    plotCompressRate();
    plotNumPoint();
    plotNumSample();
    return 0;
}
