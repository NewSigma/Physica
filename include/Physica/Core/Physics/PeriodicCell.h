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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim>
    class PeriodicCell {
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid ScalarType");
        static_assert(Dim == 2 || Dim == 3, "[Error]: Unsupported dimention");
    public:
        using LatticeMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dim, Dim>;
        using PositionMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, Dynamic, Dim>;
    protected:
        LatticeMatrix lattice;
        PositionMatrix pos;
    public:
        PeriodicCell();
        PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_);
        PeriodicCell(const PeriodicCell&) = default;
        PeriodicCell(PeriodicCell&&) noexcept = default;
        ~PeriodicCell() = default;
        /* Operators */
        PeriodicCell& operator=(PeriodicCell cell) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] const PositionMatrix& getPos() const noexcept { return pos; }
        /* Helper */
        void swap(PeriodicCell& cell) noexcept;
    };

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>::PeriodicCell()
            : lattice(LatticeMatrix::unitMatrix(Dim))
            , pos() {}

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>::PeriodicCell(LatticeMatrix lattice_, PositionMatrix pos_)
            : lattice(std::move(lattice_))
            , pos(std::move(pos_)) {}

    template<class ScalarType, unsigned int Dim>
    PeriodicCell<ScalarType, Dim>& PeriodicCell<ScalarType, Dim>::operator=(PeriodicCell cell) noexcept {
        swap(cell);
        return *this;
    }

    template<class ScalarType, unsigned int Dim>
    void PeriodicCell<ScalarType, Dim>::swap(PeriodicCell& cell) noexcept {
        lattice.swap(cell.lattice);
        pos.swap(cell.pos);
    }
}
