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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Tridiagonalization.h"

using namespace Physica::Core;

template<class MatrixType>
bool doTest(const MatrixType& source, double tolerance) {
    Tridiagonalization tri(source);
    MatrixType T = tri.getMatrixT();
    MatrixType Q = tri.getMatrixQ();
    MatrixType A = (Q * T).compute() * Q.transpose().conjugate();
    if (!matrixNear(A, source, tolerance))
        return false;
    return true;
}

int main() {
    using RealType = Scalar<Double, false>;
    {
        using MatrixType = DenseMatrix<RealType, DenseMatrixOption::Column | DenseMatrixOption::Vector, 4, 4>;
        const MatrixType temp{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        const MatrixType mat = temp + temp.transpose();
        if (!doTest(mat, 1E-14))
            return 1;
    }
    return 0;
}
