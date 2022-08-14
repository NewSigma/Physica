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

#include "Physica/Core/Physics/CrystalCell.h"

namespace Physica::Core {
    class MDCell final : PeriodicCell<Scalar<Float, false>, 3> {
    public:
        using Base = PeriodicCell<Scalar<Float, false>, 3>;
        using MassVector = Vector<ScalarType>;
    private:
        MassVector massVec;
    public:
        MDCell(CrystalCell cell);
        MDCell(Base cell, MassVector massVec_);
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const { return pos.getRow(); }
        [[nodiscard]] const MassVector& getMassVec() const { return massVec; }
        [[nodiscard]] ScalarType getMass(size_t particleID) const { return massVec[particleID]; }
        [[nodiscard]] constexpr static Type getType() noexcept { return Base::Type::Cartesian; }
    };
}
