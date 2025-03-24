/*
 * Copyright 2024-2025 Weibo He.
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
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Statistics/ProbDistribution.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

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
