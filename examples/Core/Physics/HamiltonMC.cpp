/*
 * Copyright 2025-2026 Weibo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <QApplication>
#include <print>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/MC/HamiltonMC.h"
#include "Physica/Gui/Plot/Plot.h"

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
    };
}

namespace Physica {
    template<>
    class Traits<ProbHamiltion> {
    public:
        constexpr static bool IsPeriodBoundary = false;
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
        auto p = Array<Vector2D<T>>::generate(N, [&](size_t i) {
            auto acc = hmc.step<RandomSource>(ProbHamiltion{});
            acceptV.toNextVariance(acceptM, i, acc);
            return hmc.getSample();
        });
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
