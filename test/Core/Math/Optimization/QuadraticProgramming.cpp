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
#include "Physica/Core/Math/Optimization/QuadraticProgramming/QuadraticProgramming.h"

using namespace Physica::Core;

int main() {
    using ScalarType = Scalar<Double, false>;
    {
        using VectorType = Vector<ScalarType, 2>;
        const DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector, 2, 2> G{{2, 0}, {0, 2}};
        const VectorType c{-2, -5};
        const DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> equality{};
        const DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector, 5, 3> inequality{{1, -2, -2}, {-1, -2, -6}, {-1, 2, -2}, {1, 0, 0}, {0, 1, 0}};
        const VectorType answer{1.4, 1.7};

        const VectorType initial{2, 0};
        QuadraticProgramming<ScalarType> QP(G, c, equality, inequality, initial);
        QP.compute();
        if (!vectorNear(QP.getSolution(), answer, 1E-15))
            return 1;
    }
    return 0;
}