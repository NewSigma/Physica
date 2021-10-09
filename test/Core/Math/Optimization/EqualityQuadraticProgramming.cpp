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
#include "TestHelper.h"
#include "Physica/Core/Math/Optimization/QuadraticProgramming/EqualityQuadraticProgramming.h"

using namespace Physica::Core;

int main() {
    using ScalarType = Scalar<Double, false>;
    const DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector, 3, 3> G{{6, 2, 1}, {2, 5, 2}, {1, 2, 4}};
    const Vector<ScalarType, 3> c{-8, -3, -3};
    const DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector, 2, 4> constraints{{1, 0, 1, 3}, {0, 1, 1, 0}};
    const Vector<ScalarType, 3> answer{2, -1, 1};

    const Vector<ScalarType, 3> initial{1, 1, 1};
    EqualityQuadraticProgramming<ScalarType> QP(G, c, constraints, initial);
    if (!vectorNear(QP.getSolution(), answer, 1E-8))
        return 1;
    return 0;
}
