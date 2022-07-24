/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/FokkerPlanckModel.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/Rectangle1.h"
#include <fstream>
#include "Physica/Core/IO/VTKFile.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType, 2>;

ScalarType force(VectorType pos) {
    const ScalarType x = pos[0];
    const ScalarType p = pos[1];
    ScalarType result = 0;
    if (std::isfinite(x.getTrivial()))
        result -= x;
    result -= p * 0.1;
    return result;
}

ScalarType diffuseD(VectorType pos) {
    return ScalarType(0.1);
}

ScalarType d_diffuseD(VectorType pos) {
    return ScalarType(0);
}

ScalarType initial(VectorType pos) {
    constexpr double square_sigma = 1.0 / 9;
    constexpr double square_sigma_2 = 1 / (square_sigma * 2);
    constexpr double factor = square_sigma_2 / M_PI;
    constexpr double x0 = 5;
    constexpr double p0 = 5;
    const ScalarType temp = exp(-(square(pos[0] - x0) + square(pos[1] - p0)) * square_sigma_2);
    return temp * factor;
}
/**
 * Reference:
 * [1] Spencer B F, Bergman L A. On the numerical solution of the Fokker-Planck equation for nonlinear stochastic systems[J]. Nonlinear Dynamics, 1993, 4(4):357-372.
 * [2] Andreas, Dechant, Shalom, et al. Heavy-tailed phase-space distributions beyond Boltzmann-Gibbs: Confined laser-cooled atoms in a nonthermal state[J]. Physical Review E, 2016, 94(2):022151.
 */
int main() {
    using ElementType = Rectangle1<ScalarType>;
    auto mesh = ElementType::rectangle({-M_PI_2, -M_PI_2}, {M_PI_2, M_PI_2}, 501, 501);
    mesh.addDirichletBoundary([](VectorType p) { return scalarNear(abs(p[0]), ScalarType(1), 1E-15)
                                                     || scalarNear(abs(p[1]), ScalarType(1), 1E-15); },
                              []([[maybe_unused]] VectorType p) { return ScalarType(0); });

    FokkerPlanckModel model(std::move(mesh), force, diffuseD, d_diffuseD, 1E-3, 1);
    model.setInitialCond(initial);

    for (int i = 0; i < 125; ++i)
        model.step();

    VTKFile vtk(model.getMesh(), "");
    std::ofstream fout("a.vtk");
    fout << vtk;
    return 0;
}
