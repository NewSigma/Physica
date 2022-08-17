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
#pragma once

#include "Physica/Core/Physics/PeriodicCell.h"

namespace Physica::Core {
    class CrystalCell;

    class MDCell : public PeriodicCell<Scalar<Float, false>, 3> {
        using ScalarType_ = Scalar<Float, false>;
    public:
        using Base = PeriodicCell<ScalarType_, 3>;
        using MassVector = Vector<ScalarType_>;
    private:
        using typename Base::InvLatticeMatrix;
        MassVector massVec;
        InvLatticeMatrix invLattice;
    public:
        MDCell(CrystalCell cell);
        MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_);
        /* Operations */
        void scale(ScalarType_ factor);
        void normalizeCell();
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const { return pos.getRow(); }
        [[nodiscard]] const MassVector& getMassVec() const { return massVec; }
        [[nodiscard]] ScalarType_ getMass(size_t particleID) const { return massVec[particleID]; }
        [[nodiscard]] constexpr static Type getType() noexcept { return Base::Type::Cartesian; }
    protected:
        void toDirect() { Base::toDirect(invLattice); }
        void toDirect(PositionMatrix& target) const { Base::toDirect(target, invLattice); }
    };
}
