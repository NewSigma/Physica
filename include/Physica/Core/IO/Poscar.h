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

#include <iosfwd>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    class Poscar {
    public:
        enum PoscarType : bool {
            Direct,
            Cartesian
        };

        enum CrystalSystem : char {
            Triclinic,
            Monoclinic,
            Orthohombic,
            Tetragonal,
            Hexagonal,
            Rhombohedral,
            Cubic
        };

        using ScalarType = Scalar<Float, false>;
        using LatticeMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, 3, 3>;
        using PositionMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, Dynamic, 3>;
    private:
        LatticeMatrix lattice;
        PositionMatrix pos;
        Utils::Array<size_t> numOfEachType;
        PoscarType type;
    public:
        Poscar();
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const Poscar& poscar);
        friend std::istream& operator>>(std::istream& is, Poscar& poscar);
        /* Operations */
        void standrizeLattice();
        void extendInZ(ScalarType factor);
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] LatticeMatrix getReciprocal() const noexcept;
        [[nodiscard]] CrystalSystem getCrystalSystem(double precision) const noexcept;
        [[nodiscard]] size_t getAtomCount() const noexcept;
    };
}
