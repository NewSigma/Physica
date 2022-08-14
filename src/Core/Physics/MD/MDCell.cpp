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
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    MDCell::MDCell(CrystalCell cell) : Base(std::move(cell)) {
        massVec.resize(getNumParticle());
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const auto atomicNum = cell.getAtomicNumber(i);
            massVec[i] = PhyConst<AU>::atomMass(atomicNum);
        }
    }
    MDCell::MDCell(Base cell, MassVector massVec_)
            : Base(std::move(cell))
            , massVec(std::move(massVec_))
}
