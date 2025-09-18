/*
 * Copyright 2022 Weibo He.
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
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/MDCell.h"

using namespace Physica;
using ScalarType = float64;
using MDCellType = MDCell<ScalarType>;
using CrystalCellType = CrystalCell<ScalarType>;
using LatticeMatrix = CrystalCellType::LatticeMatrix;
using PositionMatrix = CrystalCellType::PositionMatrix;

namespace {
    bool isMDCellNear(const MDCellType& cell1, const MDCellType& cell2, double precision) {
        if (!matrixNear(cell1.getLattice(), cell2.getLattice(), precision))
            return false;
        if (!matrixNear(cell1.getPos(), cell2.getPos(), precision))
            return false;
        return true;
    }
}

int main() {
    const LatticeMatrix lattice{1, 0, 0, 2, 3, 0, 4, 5, 6};
    const PositionMatrix pos{0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.5, 0.5, 0.5};

    const CrystalCellType cell1({lattice, pos, CrystalCellType::Type::Direct}, {1, 1, 2});
    CrystalCellType cell2 = cell1;
    cell2.toCartesian();

    MDCellType md1(cell1);
    MDCellType md2(cell2);
    if (!isMDCellNear(md1, md2, 1E-15))
        return 1;

    md1.normalize();
    if (!isMDCellNear(md1, md2, 1E-6))
        return 1;
    return 0;
}
