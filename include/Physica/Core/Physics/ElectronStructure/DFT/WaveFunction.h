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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    template<class ScalarType>
    class WaveFunction {
        using LatticeMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, 3, 3>;
        using Dim = std::tuple<ssize_t, ssize_t, ssize_t>;

        LatticeMatrix reciprocalLattice;
        /**
         * Coeff makes a grid of size [-baseDimX, baseDimX] * [-baseDimY, baseDimY] * [-baseDimZ, baseDimZ]
         */
        Utils::Vector<ScalarType> baseCoeffs;
        size_t baseDimX;
        size_t baseDimY;
        size_t baseDimZ;
    public:
        WaveFunction(ScalarType cutEnergy, LatticeMatrix mat);
        /* Getters */
        [[nodiscard]] size_t getPlainWaveCount() const noexcept { return baseCoeffs.getLength(); }
        [[nodiscard]] Vector<ScalarType, 3> getBaseFunc(size_t coeffIndex) const noexcept;
        [[nodiscard]] size_t dimToIndex(ssize_t x, ssize_t y, ssize_t z) const noexcept;
        [[nodiscard]] Dim indexToDim(size_t index) const noexcept;
    };

    template<class ScalarType>
    WaveFunction<ScalarType>::WaveFunction(ScalarType cutEnergy, LatticeMatrix mat) : reciprocalLattice(std::move(mat)) {
        const ScalarType factor = ScalarType(2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck);
        const ScalarType maxMoment = sqrt(PhyConst<AU>::evToHartree(cutEnergy) * factor);
        baseDimX = size_t((reciprocalLattice.row(0).norm() / maxMoment).getTrivial());
        baseDimY = size_t((reciprocalLattice.row(1).norm() / maxMoment).getTrivial());
        baseDimZ = size_t((reciprocalLattice.row(2).norm() / maxMoment).getTrivial());
        
        const size_t baseSize = (baseDimX * 2 + 1) * (baseDimY * 2 + 1) * (baseDimZ * 2 + 1);
        baseCoeffs.resize(baseSize);
    }

    template<class ScalarType>
    Vector<ScalarType, 3> WaveFunction<ScalarType>::getBaseFunc(size_t coeffIndex) const noexcept {
        const Dim dim = indexToDim(coeffIndex);
        auto[n1, n2, n3] = dim;
        return reciprocalLattice.row(0) * ScalarType(n1) +
               reciprocalLattice.row(1) * ScalarType(n2) +
               reciprocalLattice.row(2) * ScalarType(n3);
    }

    template<class ScalarType>
    size_t WaveFunction<ScalarType>::dimToIndex(ssize_t x, ssize_t y, ssize_t z) const noexcept {
        const size_t sizeX = 2 * baseDimX + 1;
        const size_t sizeXY = sizeX * (2 * baseDimY + 1);
        
        const size_t normalized_x = x + static_cast<ssize_t>(baseDimX);
        const size_t normalized_y = y + static_cast<ssize_t>(baseDimY);
        const size_t normalized_z = z + static_cast<ssize_t>(baseDimZ);
        return normalized_z * sizeXY + normalized_y * sizeX + normalized_x;
    }

    template<class ScalarType>
    Dim WaveFunction<ScalarType>::indexToDim(size_t index) const noexcept {
        const size_t sizeX = 2 * baseDimX + 1;
        const size_t sizeXY = sizeX * (2 * baseDimY + 1);

        const size_t normalized_z = index / sizeXY;
        index -= sizeXY * normalized_z;
        const size_t normalized_y = index / sizeX;
        index -= sizeX * normalized_y;
        const size_t normalized_x = index;

        return {static_cast<ssize_t>(normalized_x) - static_cast<ssize_t>(baseDimX),
                static_cast<ssize_t>(normalized_y) - static_cast<ssize_t>(baseDimY),
                static_cast<ssize_t>(normalized_z) - static_cast<ssize_t>(baseDimZ)}
    }
}
