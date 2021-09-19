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
#include <iostream>
#include "Physica/Core/Math/Optimization/GeneAlgorithm.h"
#include "Physica/Core/Physics/ElectronStructure/HF/RHFSolver.h"
#include "Physica/Core/Physics/ElectronStructure/HF/GTOnG.h"

using namespace Physica::Core;
using namespace Physica::Core::Physics;

using ScalarType = Scalar<Double, false>;

ScalarType scf_loop_He(const Vector<ScalarType, 6>& alphas) {
    Molecular<ScalarType> H2 = Molecular<ScalarType>(1);
    auto& atoms = H2.getAtoms();
    const Vector<ScalarType, 3> pos_H1{0, 0, 0};
    atoms[0] = pos_H1;
    auto& atomicNumbers = H2.getAtomicNumbers();
    atomicNumbers[0] = 1;

    ElectronConfig config = ElectronConfig(1);
    config.setOrbitState(0, ElectronConfig::SingleOccupacy);
    RHFSolver<GaussBase<ScalarType>> solver = RHFSolver<GaussBase<ScalarType>>(H2, config, 6);
    auto& baseSet = solver.getBaseSet();
    size_t i = 0;
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[0]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[1]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[2]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[3]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[4]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[5]), 0, 0, 0);
    if (!solver.compute(1E-5, 128))
        return ScalarType(1E5);
    return solver.getSelfConsistentEnergy();
}

int main() {
    srand(time(nullptr));
    const Vector<ScalarType, 6> lower{0, 0, 0, 0, 0, 0};
    const Vector<ScalarType, 6> upper{100, 100, 100, 100, 100, 100};
    GeneAlgorithm<ScalarType (*)(const Vector<ScalarType, 6>&), Vector<ScalarType, 6>> ga(scf_loop_He, {0.6, 0.1, 1000, 500, lower, upper});
    ga.solve();
    std::cout << ga.getOptimizedValue() << '\n';
    std::cout << ga.getOptimizedParams() << '\n';
    return 0;
}
