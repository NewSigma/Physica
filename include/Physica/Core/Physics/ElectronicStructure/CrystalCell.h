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

#include <unordered_set>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/IO/Poscar.h"

namespace Physica::Core {
    class ReciprocalCell;

    class CrystalCell {
    public:
        using ScalarType = Scalar<Float, false>;
        using LatticeMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, 3, 3>;
        using PositionMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, Dynamic, 3>;
    private:
        LatticeMatrix lattice;
        PositionMatrix pos;
        Utils::Array<uint16_t> atomicNumbers;
    public:
        CrystalCell(LatticeMatrix lattice_, PositionMatrix pos_, Utils::Array<uint16_t> atomicNumbers_);
        CrystalCell(Poscar poscar);
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] const PositionMatrix& getPos() const noexcept { return pos; }
        [[nodiscard]] size_t getAtomCount() const noexcept { return pos.getRow(); }
        [[nodiscard]] uint16_t getAtomicNumber(size_t ionIndex) const { return atomicNumbers[ionIndex]; }
        [[nodiscard]] ReciprocalCell reciprocal() const noexcept;
        [[nodiscard]] ScalarType getVolume() const noexcept;
        [[nodiscard]] std::unordered_set<uint16_t> getSpecies() const noexcept;
        [[nodiscard]] size_t getElectronCount() const;
    };
}
