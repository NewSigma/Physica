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

#include "PWBaseWave.h"

namespace Physica::Core::Internal {
    template<class ScalarType>
    class AbstractKSSolver {
    public:
        using ComplexType = ComplexScalar<ScalarType>;
        using Vector3D = Vector<ScalarType, 3>;
        using HermiteMatrix = DenseHermiteMatrix<ComplexType>;
        using KSOrbit = PWBaseWave<ScalarType>;
        using KSOrbitArray = Utils::Array<KSOrbit>;
        using MatrixType = DenseMatrix<ComplexType>;
        using UncenteredGrid = Grid3D<ScalarType, false>;
        using UnsignedDim = typename UncenteredGrid::Dim;
        using CenteredGrid = Grid3D<ComplexType, true>;
        using SignedDim = typename CenteredGrid::Dim;

        constexpr static size_t DIISBufferSize = 3;
        using DIISBuffer = Utils::Array<UncenteredGrid, DIISBufferSize - 1>;
        using DIISMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, DIISBufferSize, DIISBufferSize>;
    };
}
