/*
 * Copyright 2022-2025 Weibo He.
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
#include <fstream>
#include "Physica/Core/Math/Calculus/PDE/FEM/FokkerPlanckModel.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/Rectangle1.h"
#include "Physica/Core/IO/VTKFile.h"

using namespace Physica;
using T = float64;
constexpr double gammaY = 0.1;
constexpr double diffuseD_ = 0.1;
constexpr double massM = 1;

namespace {
    T force(Vector2D<T> pos) {
        const T x = pos[0];
        const T p = pos[1];
        return -x - p * gammaY;
    }

    T diffuseD([[maybe_unused]] Vector2D<T> pos) {
        return T(diffuseD_);
    }

    T d_diffuseD([[maybe_unused]] Vector2D<T> pos) {
        return T(0);
    }

    T initial(Vector2D<T> pos) {
        constexpr double square_sigma = 1.0 / 9;
        constexpr double square_sigma_2 = 1 / (square_sigma * 2);
        constexpr double factor = square_sigma_2 / M_PI;
        constexpr double x0 = 5;
        constexpr double p0 = 5;
        const T temp = exp(-(square(pos[0] - x0) + square(pos[1] - p0)) * square_sigma_2);
        return temp * factor;
    }

    T theory_stationary(Vector2D<T> pos) {
        const T a = T(gammaY / (diffuseD_ * 2));
        const T factor = a * std::sqrt(massM) / M_PI;
        return exp(-a * (square(pos[0]) * massM + square(pos[1]))) * factor;
    }
}
/**
 * Reference:
 * [1] Nonlinear Dyn 4, 357–372 (1993); https://doi.org/10.1007/BF00120671
 */
int main() {
    using ElementType = Rectangle1<T>;
    auto mesh = ElementType::rectangle({-10, -10}, {10, 10}, 201, 201);
    mesh.addDirichletBoundary([](Vector2D<T> p) { return scalarNear(abs(p[0]), T(10), 1E-12)
                                                     || scalarNear(abs(p[1]), T(10), 1E-12); },
                              []([[maybe_unused]] Vector2D<T> p) { return T(0); });

    FokkerPlanckModel model(std::move(mesh), force, diffuseD, d_diffuseD, 1E-3, massM);
    model.setInitialCond(initial);

    std::array<char, 16> buffer{}; //16 is enough
    for (int i = 0; i < 10000; ++i) {
        if (i % 250 == 0) {
            std::format_to_n(buffer.data(), buffer.size(), "FP_%d.vtk", i / 250);
            std::ofstream fout(buffer.data());
            fout << VTKFile<ElementType>(model.getMesh(), "FokkerPlanck");
        }
        model.step();
    }

    /* Test stationary */ {
        for (int i = 0; i < 90000; ++i)
            model.step();
        if (!scalarNear(model({0, 0}), theory_stationary({0, 0}), 1E-3))
            return 1;
    }
    return 0;
}
