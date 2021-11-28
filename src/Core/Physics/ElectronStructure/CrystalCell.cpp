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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"
#include "Physica/Core/Physics/ElectronStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronStructure/ReciprocalCell.h"

namespace Physica::Core {
    CrystalCell::CrystalCell(LatticeMatrix lattice_, PositionMatrix pos_, Utils::Array<unsigned char> atomicNumbers_)
            : lattice(std::move(lattice_))
            , pos(std::move(pos_))
            , atomicNumbers(std::move(atomicNumbers_)) {}

    CrystalCell::CrystalCell(Poscar poscar) : lattice(poscar.getLattice()), pos(poscar.getPos()) {}

    ReciprocalCell CrystalCell::reciprocal() const noexcept {
        LatticeMatrix result{};
        result.row(0) = lattice.row(1).crossProduct(lattice.row(2));
        result.row(1) = lattice.row(2).crossProduct(lattice.row(0));
        result.row(2) = lattice.row(0).crossProduct(lattice.row(1));
        const ScalarType volume = abs(lattice.row(0).asVector() * result.row(0));
        const ScalarType factor = ScalarType(2 * M_PI) / volume;
        result *= factor;
        return ReciprocalCell(std::move(result));
    }

    typename CrystalCell::ScalarType CrystalCell::getVolume() const noexcept {
        return abs((lattice.row(0).crossProduct(lattice.row(1))).compute() * lattice.row(2));
    }
}
