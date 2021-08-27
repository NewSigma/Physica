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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/RealSchur.h"

using namespace Physica::Core;

int main() {
    {
        using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3>;
        const MatrixType mat{{-149, 537, -27}, {-50, 180, -9}, {-154, 546, -25}};
        const MatrixType answer{{1, 0, 0}, {7.111887749987, 2, 0}, {815.8705908737, -55.02363128693, 3}};
        RealSchur schur(mat, true);
        if (!matrixNear(schur.getMatrixT(), answer, 1E-11))
            return 1;
        MatrixType A = schur.getMatrixU() * (schur.getMatrixT() * schur.getMatrixU().transpose()).compute();
        if (!matrixNear(A, mat, 1E-15))
            return 1;
    }
    return 0;
}
