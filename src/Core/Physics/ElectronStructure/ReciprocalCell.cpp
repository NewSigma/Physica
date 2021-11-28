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
#include "Physica/Core/Physics/ElectronStructure/ReciprocalCell.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"

namespace Physica::Core {
    ReciprocalCell::ReciprocalCell(LatticeMatrix lattice_) : lattice(std::move(lattice_)) {}

    typename ReciprocalCell::ScalarType ReciprocalCell::getMinNorm() const noexcept {
        const ScalarType norm1 = lattice.row(0).norm();
        const ScalarType norm2 = lattice.row(1).norm();
        const ScalarType norm3 = lattice.row(2).norm();
        return norm1 > norm2 ? (norm2 > norm3 ? norm3 : norm2)
                             : (norm1 > norm3 ? norm3 : norm1);
    }

    typename ReciprocalCell::ScalarType ReciprocalCell::getVolume() const noexcept {
        return abs((lattice.row(0).crossProduct(lattice.row(1))).compute() * lattice.row(2));
    }
}
