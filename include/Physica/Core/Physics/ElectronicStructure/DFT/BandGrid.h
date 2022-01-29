/*
 * Copyright 2021-2022 WeiBo He.
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
#include "KPointInfo.h"
#include "Grid3D.h"

namespace Physica::Core {
    template<class ScalarType, bool isSpinPolarized> class BandGrid;

    template<class ScalarType, bool isSpinPolarized>
    class BandGrid : public Grid3D<KPointInfo<ScalarType, isSpinPolarized>, true> {
        using Base = Grid3D<KPointInfo<ScalarType, isSpinPolarized>, true>;
        using typename Base::LatticeMatrix;
    public:
        using KPoint = Vector<ScalarType, 3>;
    public:
        template<class MatrixType>
        BandGrid(ScalarType cutEnergy, const LValueMatrix<MatrixType>& repLatt, size_t dimX, size_t dimY, size_t dimZ);
    private:
        template<class MatrixType>
        static LatticeMatrix getGridLattice(const LValueMatrix<MatrixType>& repLatt, size_t dimX, size_t dimY, size_t dimZ);
    };

    template<class ScalarType, bool isSpinPolarized>
    template<class MatrixType>
    BandGrid<ScalarType, isSpinPolarized>::BandGrid(ScalarType cutEnergy,
                                                    const LValueMatrix<MatrixType>& repLatt,
                                                    size_t dimX,
                                                    size_t dimY,
                                                    size_t dimZ)
            : Base(getGridLattice(repLatt, dimX, dimY, dimZ), dimX, dimY, dimZ, Base::sizeFromCutEnergy(cutEnergy, repLatt)) {}

    template<class ScalarType, bool isSpinPolarized>
    template<class MatrixType>
    typename BandGrid<ScalarType, isSpinPolarized>::LatticeMatrix
    BandGrid<ScalarType, isSpinPolarized>::getGridLattice(const LValueMatrix<MatrixType>& repLatt, size_t dimX, size_t dimY, size_t dimZ) {
        using LatticeScalar = typename LatticeMatrix::ScalarType;
        LatticeMatrix result{};
        result.row(0) = repLatt.row(0).asVector() * reciprocal(LatticeScalar((dimX + 1) * 2));
        result.row(1) = repLatt.row(1).asVector() * reciprocal(LatticeScalar((dimY + 1) * 2));
        result.row(2) = repLatt.row(2).asVector() * reciprocal(LatticeScalar((dimZ + 1) * 2));
        return result;
    }
}
