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
#include "TestHelper.h"
#include "Physica/Core/Physics/ElectronStructure/HF/RHFSolver.h"
#include "Physica/Core/Physics/ElectronStructure/HF/GTOnG.h"

using namespace Physica::Core;
using namespace Physica::Core::Physics;

using ScalarType = Scalar<Double, false>;

ScalarType scf_solve(size_t atomicNumber, const ElectronConfig& config, const Vector<ScalarType>& alphas) {
    Molecular<ScalarType> system = Molecular<ScalarType>(1);
    auto& atoms = system.getAtoms();
    const Vector<ScalarType, 3> pos_atom{0, 0, 0};
    atoms[0] = pos_atom;
    auto& atomicNumbers = system.getAtomicNumbers();
    atomicNumbers[0] = atomicNumber;

    RHFSolver<GaussBase<ScalarType>> solver = RHFSolver<GaussBase<ScalarType>>(system, config, alphas.getLength());
    auto& baseSet = solver.getBaseSet();
    for (size_t i = 0; i < alphas.getLength(); ++i)
        baseSet[i] = GaussBase<ScalarType>(pos_atom, abs(alphas[i]), 0, 0, 0);
    if (!solver.compute(1E-5, 128))
        return ScalarType(1E5);
    return solver.getSelfConsistentEnergy();
}
/**
 * Reference:
 * [1] Koga T, Tatewaki H, Shimazaki T. Chemically reliable uncontracted Gaussian-type basis sets for atoms H to Lr [J].Chemical Physics Letters, 2000, 328(4-6):473-482.
 */
int main() {
    {
        ElectronConfig config = ElectronConfig(1);
        config.setOrbitState(0, ElectronConfig::SingleOccupacy);
        const Vector<ScalarType, 6> base{8.2921890E1, 1.2452437E1, 2.8330562, 8.0001038E-1, 2.5859469E-1, 8.9968966E-2};
        if (!scalarNear(scf_solve(1, config, base), ScalarType(-0.499945570), 1E-9))
            return 1;
    }
    return 0;
}
