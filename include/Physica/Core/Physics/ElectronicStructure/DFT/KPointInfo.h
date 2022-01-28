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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"

namespace Physica::Core {
    template<class ScalarType>
    class KPointInfo {
        using ComplexType = ComplexScalar<ScalarType>;
        using RawDataType = EigenSolver<DenseMatrix<ComplexType>>;

        RawDataType bandUp;
        RawDataType bandDown;
    public:
        KPointInfo(size_t plainWaveCount);
        /* Setters */
        void setRawData(RawDataType& newBandUp, RawDataType& newBandDown);
        /* Getters */
        const RawDataType& getBandUp() const noexcept { return bandUp; }
        const RawDataType& getBandDown() const noexcept { return bandDown; }
    };

    template<class ScalarType>
    KPointInfo<ScalarType>::KPointInfo(size_t plainWaveCount) : bandUp(plainWaveCount), bandDown(plainWaveCount) {}

    template<class ScalarType>
    void KPointInfo<ScalarType>::setRawData(RawDataType& newBandUp, RawDataType& newBandDown) {
        std::swap(bandUp, newBandUp);
        std::swap(bandDown, newBandDown);
    }
}
