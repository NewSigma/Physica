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

struct InitialIntegrator {
    template<class Functor>
    static ScalarType run(Functor func) {
        return func({0, 0});
    }
};

struct Integrator {
    template<class Functor>
    static ScalarType run(Functor func) {
        Integrate<Rectangular, ScalarType, 1> integrator(IntegrateRange<ScalarType, 1>({-1}, {1}), 0.2);
        return integrator.solve([&](ScalarType x) {
            return integrator.solve([=](ScalarType y) -> ScalarType { return func({x, y}); });
        });
    }
};

ScalarType force(VectorType pos) {
    constexpr double gammaY = 1;
    const ScalarType temp1 = -square(pos[0]) + 1;
    const ScalarType temp2 = -square(pos[1]) + 1;
    if (temp1 < std::numeric_limits<ScalarType>::epsilon() || temp2 < std::numeric_limits<ScalarType>::epsilon())
        return ScalarType::Zero();
    const ScalarType x = pos[0] / sqrt(temp1);
    const ScalarType p = pos[1] / sqrt(temp2);
    return -x / (ScalarType(1) + square(x)) - p / (ScalarType(1) + square(p)) * gammaY;
}

ScalarType diffuseD(VectorType pos) {
    return reciprocal(square(pos[1]) + 1) + 1;
}

ScalarType d_diffuseD(VectorType pos) {
    return ScalarType(-2) * pos[1] / (square(square(pos[1]) + 1));
}

ScalarType initial(VectorType pos) {
    if (scalarNear(pos.norm(), ScalarType(0), 1E-15))
        return ScalarType(1);
    return ScalarType(0);
}
/**
 * Reference:
 * [1]Andreas, Dechant, Shalom, et al. Heavy-tailed phase-space distributions beyond Boltzmann-Gibbs: Confined laser-cooled atoms in a nonthermal state[J]. Physical Review E, 2016, 94(2):022151.
 */
int main() {
    using ElementType = Rectangle1<ScalarType>;
    auto mesh = ElementType::rectangle({-1, -1}, {1, 1}, 41, 41);
    mesh.addDirichletBoundary([](VectorType p) { return scalarNear(abs(p[0]), ScalarType(1), 1E-15)
                                                     || scalarNear(abs(p[1]), ScalarType(1), 1E-15); },
                              []([[maybe_unused]] VectorType p) { return ScalarType(0); });

    FokkerPlanckModel model(std::move(mesh), force, diffuseD, d_diffuseD, 0.01, 1);
    model.setInitialCond<InitialIntegrator>(initial);

    for (int i = 0; i < 10; ++i)
        model.step<Integrator>();

    VTKFile vtk(model.getMesh(), "");
    std::ofstream fout("a.vtk");
    fout << vtk;
    return 0;
}
