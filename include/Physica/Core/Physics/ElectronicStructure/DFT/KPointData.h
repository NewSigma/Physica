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
    class KPointData {
        using ComplexType = ComplexScalar<ScalarType>;
        using DataType = EigenSolver<DenseMatrix<ComplexType>>;

        DataType bandUp;
        DataType bandDown;
    public:
        KPointData(size_t plainWaveCount);
        /* Setters */
        void setData(DataType& newBandUp, DataType& newBandDown);
    };

    template<class ScalarType>
    KPointData<ScalarType>::KPointData(size_t plainWaveCount) : bandUp(plainWaveCount), bandDown(plainWaveCount) {}

    template<class ScalarType>
    void KPointData<ScalarType>::setData(DataType& newBandUp, DataType& newBandDown) {
        std::swap(bandUp, newBandUp);
        std::swap(bandDown, newBandDown);
    }
}
