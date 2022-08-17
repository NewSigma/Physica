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
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"
#include <iostream>

namespace Physica::Core {
    MDCell::MDCell(CrystalCell cell)
            : Base(std::move(cell))
            , invLattice(Base::makeInvLattice()) {
        massVec.resize(getNumParticle());
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const auto atomicNum = cell.getAtomicNumber(i);
            massVec[i] = PhyConst<AU>::atomMass(atomicNum);
        }
    }

    MDCell::MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_)
            : Base(std::move(lattice), std::move(pos))
            , massVec(std::move(massVec_))
            , invLattice(Base::makeInvLattice()) {}

    void MDCell::scale(ScalarType_ factor) {
        assert(factor.isPositive());
        lattice *= factor;
        pos *= factor;
        invLattice *= Core::reciprocal(factor);
    }

    void MDCell::normalizeCell() {
        Base::toDirect(invLattice);
        for (auto& elem : pos) {
            const int integer = float(elem);
            elem -= ScalarType_(integer - elem.isNegative());
            assert(ScalarType_::Zero() <= elem && elem <= ScalarType_::One());
        }
        Base::toCartesian();
    }
}
