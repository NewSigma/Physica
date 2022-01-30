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
    class BandGrid {
        using KPoint = KPointInfo<ScalarType, isSpinPolarized>;
        using KPoints = Utils::Array<KPoint>;

        Utils::Array<KPoint> kPoints;
    public:
        template<class MatrixType>
        BandGrid(ScalarType cutEnergy, const LValueMatrix<MatrixType>& repLatt, size_t kPointX, size_t kPointY, size_t kPointZ);
        /* Getters */
        [[nodiscard]] KPoints& getKPoints() { return kPoints; }
        [[nodiscard]] const KPoints& getKPoints() const noexcept { return kPoints; }
    };

    template<class ScalarType, bool isSpinPolarized>
    template<class MatrixType>
    BandGrid<ScalarType, isSpinPolarized>::BandGrid(ScalarType cutEnergy,
                                                    const LValueMatrix<MatrixType>& repLatt,
                                                    size_t kPointX,
                                                    size_t kPointY,
                                                    size_t kPointZ)
            : kPoints(kPointX * kPointY * kPointZ) {
        assert(kPoints.getLength() != 0);
        size_t kPointID = 0;
        const size_t plainWaveCount = Grid3D<double, true>::sizeFromCutEnergy(cutEnergy, repLatt); //TODO: signed/unsigned character can be moved to father class
        
        const ScalarType stepX = reciprocal(ScalarType(kPointX));
        const ScalarType stepY = reciprocal(ScalarType(kPointY));
        const ScalarType stepZ = reciprocal(ScalarType(kPointZ));

        Vector<ScalarType, 3> k{};
        ScalarType& kx = k[0];
        ScalarType& ky = k[1];
        ScalarType& kz = k[2];

        kx = (ScalarType(1) - ScalarType(kPointX)) / ScalarType(2 * kPointX);
        for (size_t x = 1; x <= kPointX; ++x) {
            ky = (ScalarType(1) - ScalarType(kPointY)) / ScalarType(2 * kPointY);
            for (size_t y = 1; y <= kPointY; ++y) {
                kz = (ScalarType(1) - ScalarType(kPointZ)) / ScalarType(2 * kPointZ);
                for (size_t z = 1; z <= kPointZ; ++z) {
                    kPoints[kPointID] = KPointInfo<ScalarType, isSpinPolarized>(k, plainWaveCount);
                    kz += stepZ;
                    ++kPointID;
                }
                ky += stepY;
            }
            kx += stepX;
        }
    }
}
