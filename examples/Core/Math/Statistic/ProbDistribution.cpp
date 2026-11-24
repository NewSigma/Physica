/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Statistics/ProbDistribution.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using RandomSource = Random<>;

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Plot* plot = new Plot(-5, 5, 0, 1, 2, 0.5);
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%.1f");
    plot->getAxisX()->setTitleText("y");
    plot->getAxisY()->setTitleText("P(y)");

    constexpr double a = 1;
    constexpr int NumSample = 100000;
    ProbDistribution<float64> dist(-5, 5, 1000);
    for (int i = 0; i < NumSample; ++i) {
        dist.sample([](int n) -> float64 {
            const auto temp = RandomSource::random_int(n, 0, 1);
            float64 result = 0;
            for (int i : temp)
                result += float64(i == 0 ? 1.0 : -1.0);
            return result * ln1p(sqrt(float64(1) - exp(-float64(a / n))));
        }(1000));
    }
    const auto x = dist.makePosition();
    auto y = dist.makeDistribution();
    y *= reciprocal(y.max());
    const VectorND<float64> y1 = exp(square(x) * float64(-0.5));
    plot->line(x, y);
    plot->line(x, y1);
    plot->show();
    return QApplication::exec();
}
