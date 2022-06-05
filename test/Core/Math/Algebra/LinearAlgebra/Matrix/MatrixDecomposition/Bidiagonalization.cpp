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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Bidiagonalization.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector>;

template<class MatrixType>
bool doTest(const MatrixType& source, double tolerance) {
    Bidiagonalization obj(source);
    MatrixType U = obj.getMatrixU();
    MatrixType V = obj.getMatrixV();
    MatrixType B = obj.getMatrixB();
    MatrixType A = (U * B).compute() * V.transpose();
    if (!matrixNear(A, source, tolerance))
        return false;
    return true;
}

int main() {
    {
        const MatrixType mat{{1, 2, 3, 4, 5}, {5, 6, 7, 8, 9}, {9, 10, 11, 12, 13}, {7, 6, -8, -9, 5}};
        if (!doTest(mat, 1E-15))
            return 1;
    }
    return 0;
}
