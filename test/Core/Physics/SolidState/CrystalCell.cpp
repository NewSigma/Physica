/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float32;
using CrystalCellType = CrystalCell<ScalarType>;
using LatticeMatrix = CrystalCellType::LatticeMatrix;
using PositionMatrix = CrystalCellType::PositionMatrix;

int main() {
    const LatticeMatrix lattice{1, 0, 0, 2, 3, 0, 4, 5, 6};
    const PositionMatrix pos{0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.5, 0.5, 0.5};
    const CrystalCellType cell_direct({lattice, pos, CrystalCellType::Type::Direct}, {1, 1, 8});

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
    const Array<uint16_t> answer_atomic{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8};
    CrystalCellType supercell_direct = cell_direct.makeSuperCell(3, 2, 1);
    expect(matrixNear(answer_lattice, supercell_direct.getLattice(), 1E-15));
    expect(matrixNear(answer_pos, supercell_direct.getPos(), 1E-15));

    for (size_t i = 0; i < answer_atomic.getLength(); ++i)
        expect(answer_atomic[i] == supercell_direct.getAtomicNumber(i));

    CrystalCellType supercell_cartesian = cell_direct.makeSuperCell(3, 2, 1);
    supercell_cartesian.toCartesian();
    expect(matrixNear(answer_lattice, supercell_cartesian.getLattice(), 1E-15));

    supercell_direct.toCartesian();
    expect(matrixNear(supercell_direct.getPos(), supercell_cartesian.getPos(), 1E-15));

    for (size_t i = 0; i < answer_atomic.getLength(); ++i)
        expect(answer_atomic[i] == supercell_cartesian.getAtomicNumber(i));
    return 0;
}
