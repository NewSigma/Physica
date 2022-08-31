/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"

namespace Physica::Core {
    template<class ScalarType>
    class ReciprocalCell {
        using LatticeMatrix = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, 3, 3>;

        LatticeMatrix lattice;
    public:
        ReciprocalCell() = default;
        ReciprocalCell(LatticeMatrix lattice_);
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept;
        [[nodiscard]] ScalarType getMinNorm() const noexcept;
        [[nodiscard]] ScalarType getVolume() const noexcept;
    };

    template<class ScalarType>
    ReciprocalCell<ScalarType>::ReciprocalCell(LatticeMatrix lattice_) : lattice(std::move(lattice_)) {}

    template<class ScalarType>
    const typename ReciprocalCell<ScalarType>::LatticeMatrix& ReciprocalCell<ScalarType>::getLattice() const noexcept {
        return lattice;
    }

    template<class ScalarType>
    ScalarType ReciprocalCell<ScalarType>::getMinNorm() const noexcept {
        const ScalarType norm1 = lattice.row(0).norm();
        const ScalarType norm2 = lattice.row(1).norm();
        const ScalarType norm3 = lattice.row(2).norm();
        return norm1 > norm2 ? (norm2 > norm3 ? norm3 : norm2)
                             : (norm1 > norm3 ? norm3 : norm1);
    }

    template<class ScalarType>
    ScalarType ReciprocalCell<ScalarType>::getVolume() const noexcept {
        return abs((lattice.row(0).crossProduct(lattice.row(1))).compute() * lattice.row(2).asVector());
    }
}
