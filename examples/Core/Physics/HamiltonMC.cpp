/*
 * Copyright 2025-2026 Weibo He.
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
#include <QApplication>
#include <print>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/MC/HamiltonMC.h"

import Physica.Gui.Plot;

using namespace Physica;
using T = float64;
using RandomSource = Random<>;

namespace {
    class ProbHamiltion : private EmptyForceModel<T, 1> {
        using Base = EmptyForceModel<T, 1>;
        using Cell = Base::MDCellType;
    public:
        template<ExecutePolicy>
        [[nodiscard]] static T potentialV(const Cell& cell) {
            const auto& pos = cell.getPos();
            T x = pos[0, 0];
            T y = pos[1, 0];
            return square(x) + T(2) * square(y) + T(2) * x * y;
        }

        template<ExecutePolicy P>
        static void forceAsync(const MDCellType& cell, Vector auto& result) {
            const auto& pos = cell.getPos();
            T x = pos[0, 0];
            T y = pos[1, 0];
            result = Vector2D<T>{T(-2) * x - T(2) * y, T(-2) * x - T(4) * y};
        }
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return false; }
    };
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Plot* plot = new Plot(-4, 4, -4, 4, 2, 2);
    plot->getLegend().show();
    plot->getLegend().setAlignment(Qt::AlignRight);
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%d");
    plot->getAxisX()->setTitleText("x");
    plot->getAxisY()->setTitleText("y");

    constexpr int N = 2000;
    VectorND<T> x(N), y(N);
    {
        auto hmc = HamiltonMC<T>({1, 1});
        hmc.warmup<RandomSource>(1000, ProbHamiltion{});

        T acceptM, acceptV;
        auto p = Array<Vector2D<T>>::generate([&](size_t i) {
            auto acc = hmc.step<RandomSource>(ProbHamiltion{});
            acceptV.toNextVariance(acceptM, i, acc);
            return hmc.getSample();
        }, N);
        std::println("{}({})", acceptM, sqrt(acceptV));

        std::ranges::transform(p, x.begin(), [](Vector2D<T> in) { return in[0]; });
        std::ranges::transform(p, y.begin(), [](Vector2D<T> in) { return in[1]; });
    }
    auto& s = plot->scatter(x, y);
    s.setBorderColor(s.color());
    s.setColor(Qt::transparent);
    s.setMarkerSize(8);

    plot->show();
    return QApplication::exec();
}
