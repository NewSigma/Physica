/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Physics/ElectronicStructure/HF/RHFSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/HF/GaussBase.h"

using namespace Physica;
using ScalarType = float64;

constexpr double precitionGoal = 1E-8;

namespace {
    ScalarType scf_solve(size_t atomicNumber, const ElectronConfig& config, const VectorND<ScalarType>& alphas, size_t maxIte) {
        Molecular<ScalarType> system = Molecular<ScalarType>(1);
        auto& atoms = system.getAtoms();
        const Vector3D<ScalarType> pos_atom{0, 0, 0};
        atoms[0] = pos_atom;
        auto& atomicNumbers = system.getAtomicNumbers();
        atomicNumbers[0] = atomicNumber;

        RHFSolver<GaussBase<ScalarType>> solver = RHFSolver<GaussBase<ScalarType>>(system, config, alphas.getLength());
        auto& baseSet = solver.getBaseSet();
        for (size_t i = 0; i < alphas.getLength(); ++i)
            baseSet[i] = GaussBase<ScalarType>(pos_atom, abs(alphas[i]), 0, 0, 0);
        if (!solver.compute(precitionGoal, maxIte))
            return ScalarType(1E5);
        return solver.getSelfConsistentEnergy();
    }
}
/**
 * Reference:
 * [1] Chemical Physics Letters 328(4-6), 473-482 (2000); https://doi.org/10.1016/S0009-2614(00)00948-9
 */
int main() {
    {
        ElectronConfig config = ElectronConfig(1);
        config.setOrbitState(0, ElectronConfig::SingleOccupacy);
        const VectorND<ScalarType> base{8.2921890E1, 1.2452437E1, 2.8330562, 8.0001038E-1, 2.5859469E-1, 8.9968966E-2};
        if (!scalarNear(scf_solve(1, config, base, 1), ScalarType(-0.499945570), precitionGoal))
            return 1;
    }
    {
        ElectronConfig config = ElectronConfig(1);
        config.setOrbitState(0, ElectronConfig::DoubleOccupacy);
        const VectorND<ScalarType> base{2.3406425E2, 3.5174075E1, 7.9911108, 2.2124201, 6.6706872E-1, 2.0894727E-1};
        if (!scalarNear(scf_solve(2, config, base, 8), ScalarType( -2.86115334), precitionGoal))
            return 1;
    }
    return 0;
}
