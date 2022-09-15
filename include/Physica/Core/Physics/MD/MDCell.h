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

    class MDCell : public PeriodicCell<Scalar<Double, false>, 3> {
    public:
        using ScalarType = Scalar<Double, false>;
        using Base = PeriodicCell<ScalarType, 3>;
        using MassVector = Vector<ScalarType>;
    private:
        using typename Base::InvLatticeMatrix;
        MassVector massVec;
        InvLatticeMatrix invLattice;
    public:
        MDCell(CrystalCell cell);
        MDCell(LatticeMatrix lattice, PositionMatrix pos, MassVector massVec_);
        /* Operations */
        void scale(ScalarType factor);
        void normalizeCell();
        void toDirect(PositionMatrix& target) const { Base::toDirect(target, invLattice); }
        /* Getters */
        [[nodiscard]] const MassVector& getMassVec() const { return massVec; }
        [[nodiscard]] ScalarType getMass(size_t particleID) const { return massVec[particleID]; }
        [[nodiscard]] constexpr static Type getType() noexcept { return Base::Type::Cartesian; }
    protected:
        void toDirect() { Base::toDirect(invLattice); }
    };
}
