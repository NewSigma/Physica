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
#include "Physica/Core/Math/Optimization/GeneAlgorithm.h"
#include "Physica/Core/Physics/ElectronStructure/HF/HFSolver.h"
#include "Physica/Core/Physics/ElectronStructure/HF/GTOnG.h"

using namespace Physica::Core;
using namespace Physica::Core::Physics;

using ScalarType = Scalar<Double, false>;

ScalarType scf_loop_He(const Vector<ScalarType, 4>& alphas) {
    Molecular<ScalarType> He = Molecular<ScalarType>(1);
    auto& atoms = He.getAtoms();
    const Vector<ScalarType> pos_He{0, 0, 0};
    atoms[0] = pos_He;
    auto& atomicNumbers = He.getAtomicNumbers();
    atomicNumbers[0] = 2;

    HFSolver<GaussBase<ScalarType>> solver = HFSolver<GaussBase<ScalarType>>(He, 2, 4);
    auto& baseSet = solver.getBaseSet();
    size_t i = 0;
    baseSet[i++] = GaussBase<ScalarType>(pos_He, abs(alphas[0]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_He, abs(alphas[1]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_He, abs(alphas[2]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_He, abs(alphas[3]), 0, 0, 0);
    if (!solver.compute(1E-5, 128))
        return ScalarType(1E5);
    return solver.getSelfConsistentEnergy();
}

int main() {
    srand(time(nullptr));
    Vector<ScalarType, 4> x{randomScalar<ScalarType>(ScalarType(0), ScalarType(10)),
                            randomScalar<ScalarType>(ScalarType(0), ScalarType(10)),
                            randomScalar<ScalarType>(ScalarType(0), ScalarType(10)),
                            randomScalar<ScalarType>(ScalarType(0), ScalarType(10))};
    const Vector<ScalarType, 4> lower{0, 0, 0, 0};
    const Vector<ScalarType, 4> upper{50, 50, 50, 50};
    GeneAlgorithm<ScalarType (*)(const Vector<ScalarType, 4>&), Vector<ScalarType, 4>> ga(scf_loop_He, {0.6, 0.1, 100, 50, lower, upper});
    ga.solve();
    std::cout << ga.getOptimizedValue() << '\n';
    std::cout << ga.getOptimizedParams() << '\n';
    return 0;
}
