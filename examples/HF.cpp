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
#include "Physica/Core/Physics/ElectronStructure/HF/HFSolver.h"
#include "Physica/Core/Physics/ElectronStructure/HF/GTOnG.h"

using namespace Physica::Core;
using namespace Physica::Core::Physics;

using ScalarType = Scalar<Double, false>;

int main() {
    {
        Molecular<ScalarType> He = Molecular<ScalarType>(1);
        auto& atoms = He.getAtoms();
        const Vector<ScalarType> pos_He{0, 0, 0};
        atoms[0] = pos_He;
        auto& atomicNumbers = He.getAtomicNumbers();
        atomicNumbers[0] = 2;

        HFSolver<GaussBase<ScalarType>> solver = HFSolver<GaussBase<ScalarType>>(He, 2, 4);
        auto& baseSet = solver.getBaseSet();
        size_t i = 0;
        baseSet[i++] = GaussBase<ScalarType>(pos_He, 0.298073, 0, 0, 0);
        baseSet[i++] = GaussBase<ScalarType>(pos_He, 1.242567, 0, 0, 0);
        baseSet[i++] = GaussBase<ScalarType>(pos_He, 5.782948, 0, 0, 0);
        baseSet[i++] = GaussBase<ScalarType>(pos_He, 38.474970, 0, 0, 0);
        if (!solver.compute(1E-5, 20)) {
            std::cout << "[Warning]: Not converged\n";
            return 1;
        }
        std::cout << "Self consistent energy(in hartree): " << solver.getSelfConsistentEnergy() << '\n';
        std::cout << "Totle energy(in hartree): " << solver.getTotalEnergy() << '\n';
    }
    {
        Molecular<ScalarType> LiF = Molecular<ScalarType>(2);
        auto& atoms = LiF.getAtoms();
        const Vector<ScalarType> pos_Li{0, 0, 0};
        const Vector<ScalarType> pos_F{0, 0, 1};
        atoms[0] = pos_Li;
        atoms[1] = pos_F;
        auto& atomicNumbers = LiF.getAtomicNumbers();
        atomicNumbers[0] = 3;
        atomicNumbers[1] = 9;

        HFSolver<GTO3G<ScalarType>> solver = HFSolver<GTO3G<ScalarType>>(LiF, 12, 12);
        auto& baseSet = solver.getBaseSet();
        size_t i = 0;
        for (; i < 3; ++i)
            baseSet[i] = GTO3G<ScalarType>::randomBase(pos_Li);
        for (; i < 7; ++i)
            baseSet[i] = GTO3G<ScalarType>::randomBase(pos_F);
        baseSet[i++] = GTO3G<ScalarType>::randomBase(pos_F, 1, 0, 0);
        baseSet[i++] = GTO3G<ScalarType>::randomBase(pos_F, 1, 0, 0);
        baseSet[i++] = GTO3G<ScalarType>::randomBase(pos_F, 0, 1, 0);
        baseSet[i++] = GTO3G<ScalarType>::randomBase(pos_F, 0, 1, 0);
        baseSet[i++] = GTO3G<ScalarType>::randomBase(pos_F, 0, 0, 1);
        if (!solver.compute(1E-5, 50000)) {
            std::cout << "[Warning]: Not converged\n";
            return 1;
        }
        std::cout << "Self consistent energy(in hartree): " << solver.getSelfConsistentEnergy() << '\n';
        std::cout << "Totle energy(in hartree): " << solver.getTotalEnergy() << '\n';
    }
    return 0;
}
