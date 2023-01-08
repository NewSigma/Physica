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
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"

using namespace Physica::Core;

int main() {
    using LatticeMatrix = typename CrystalCell::LatticeMatrix;
    using PositionMatrix = typename CrystalCell::PositionMatrix;
    const LatticeMatrix lattice{1, 0, 0, 2, 3, 0, 4, 5, 6};
    const PositionMatrix pos{0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.5, 0.5, 0.5};
    const CrystalCell cell_direct({lattice, pos, CrystalCell::Type::Direct}, {1, 1, 8});

    const LatticeMatrix answer_lattice{3, 0, 0, 4, 6, 0, 4, 5, 6};
    const PositionMatrix answer_pos{
        0.0833333333333333, 0.1250000000000000, 0.2500000000000000,
        0.0833333333333333, 0.6250000000000000, 0.2500000000000000,
        0.4166666666666666, 0.1250000000000000, 0.2500000000000000,
        0.4166666666666666, 0.6250000000000000, 0.2500000000000000,
        0.7500000000000000, 0.1250000000000000, 0.2500000000000000,
        0.7500000000000000, 0.6250000000000000, 0.2500000000000000,
        0.0833333333333333, 0.3750000000000000, 0.7500000000000000,
        0.0833333333333333, 0.8750000000000000, 0.7500000000000000,
        0.4166666666666666, 0.3750000000000000, 0.7500000000000000,
        0.4166666666666666, 0.8750000000000000, 0.7500000000000000,
        0.7500000000000000, 0.3750000000000000, 0.7500000000000000,
        0.7500000000000000, 0.8750000000000000, 0.7500000000000000,
        0.1666666666666667, 0.2500000000000000, 0.5000000000000000,
        0.1666666666666667, 0.7500000000000000, 0.5000000000000000,
        0.5000000000000000, 0.2500000000000000, 0.5000000000000000,
        0.5000000000000000, 0.7500000000000000, 0.5000000000000000,
        0.8333333333333333, 0.2500000000000000, 0.5000000000000000,
        0.8333333333333333, 0.7500000000000000, 0.5000000000000000
    };
    const Physica::Utils::Array<uint16_t> answer_atomic{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8};
    CrystalCell supercell_direct = cell_direct;
    supercell_direct.unitToSuper(3, 2, 1);
    if (!matrixNear(answer_lattice, supercell_direct.getLattice(), 1E-15))
        return 1;
    if (!matrixNear(answer_pos, supercell_direct.getPos(), 1E-15))
        return 1;
    for (size_t i = 0; i < answer_atomic.getLength(); ++i)
        if (answer_atomic[i] != supercell_direct.getAtomicNumber(i))
            return 1;

    CrystalCell supercell_cartesian = cell_direct;
    supercell_cartesian.toCartesian();
    supercell_cartesian.unitToSuper(3, 2, 1);
    if (!matrixNear(answer_lattice, supercell_cartesian.getLattice(), 1E-15))
        return 1;
    supercell_direct.toCartesian();
    if (!matrixNear(supercell_direct.getPos(), supercell_cartesian.getPos(), 1E-15))
        return 1;
    for (size_t i = 0; i < answer_atomic.getLength(); ++i)
        if (answer_atomic[i] != supercell_cartesian.getAtomicNumber(i))
            return 1;
    return 0;
}
