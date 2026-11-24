/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Physics/SolidState/PeriodicCell.h"

namespace Physica {
    /**
     * Angular function introduced in [1]
     *
     * Reference:
     * [1] Phys. Rev. Lett. 98, 146401 (2007); https://doi.org/10.1103/PhysRevLett.98.146401
     */
    template<Scalar T, bool IsSmallCell>
    class GaussAngular {
        using VectorType = VectorND<T>;
        using CellType = PeriodicCell<T, 3>;

        T paramEta;
        T distR;
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

    template<Scalar T, bool IsSmallCell>
    GaussAngular<T, IsSmallCell>::VectorType GaussAngular<T, IsSmallCell>::calc(const CellType& cell) const {
        VectorType result(cell.getNumParticle());
        if constexpr (IsSmallCell) {

        }
        else {

        }
    }
}
