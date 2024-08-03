/*
 * Copyright 2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"

namespace Physica::Core {
    /**
     * Angular function introduced in [1]
     * 
     * Reference:
     * [1] Phys. Rev. Lett. 98, 146401 (2007); https://doi.org/10.1103/PhysRevLett.98.146401
     */
    template<class ScalarType, bool IsSmallCell>
    class GaussAngular {
        using VectorType = Vector<ScalarType>;
        using CellType = PeriodicCell<ScalarType, 3>;

        ScalarType paramEta;
        ScalarType distR;
    public:
        GaussAngular();
        GaussAngular(const GaussAngular&) = default;
        GaussAngular(GaussAngular&&) noexcept = default;
        ~GaussAngular() = default;
        /* Operators */
        GaussAngular& operator=(GaussAngular obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] VectorType calc(const CellType& cell) const;
        void swap(GaussAngular& __restrict obj) noexcept;
    };

    template<class ScalarType, bool IsSmallCell>
    typename GaussAngular<ScalarType, IsSmallCell>::VectorType GaussAngular<ScalarType, IsSmallCell>::calc(const CellType& cell) const {
        VectorType result(cell.getNumParticle());
        if constexpr (IsSmallCell) {

        }
        else {

        }
    }
}
