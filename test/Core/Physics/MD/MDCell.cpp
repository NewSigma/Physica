/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    bool isMDCellNear(const MDCell<T>& cell1, const MDCell<T>& cell2, double precision) {
        if (!matrixNear(cell1.getLattice(), cell2.getLattice(), precision))
            return false;
        if (!matrixNear(cell1.getPos(), cell2.getPos(), precision))
            return false;
        return true;
    }

    void fromCrystal() {
        using CrystalCellType = CrystalCell<T>;
        using LatticeMatrix = CrystalCellType::LatticeMatrix;
        const LatticeMatrix lattice{1, 0, 0, 2, 3, 0, 4, 5, 6};

        using PositionMatrix = CrystalCellType::PositionMatrix;
        const PositionMatrix pos{0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.5, 0.5, 0.5};

        const CrystalCellType cell1({lattice, pos, CrystalCellType::Type::Direct}, {1, 1, 2});
        CrystalCellType cell2 = cell1;
        cell2.toCartesian();

        MDCell<T> md1(cell1);
        MDCell<T> md2(cell2);
        expect(isMDCellNear(md1, md2, 1E-15));

        md1.normalize();
        expect(isMDCellNear(md1, md2, 1E-6));
    }
}

int main() {
    fromCrystal();
    return 0;
}
