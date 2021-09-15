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

template<class MatrixType>
bool isUpperQuasiTriangle(const LValueMatrix<MatrixType>& m) {
    if (m.getRow() != m.getColumn())
        return false;
    for (size_t i = 0; i < m.getRow() - 1; ++i) {
        if (m(i + 1, i).isZero()) {
            auto col = m.col(i);
            if (!col.tail(i + 1).isZero())
                return false;
        }
        else if(i < m.getRow() - 2) {
            auto col1 = m.col(i);
            auto col2 = m.col(i + 1);
            if (!col1.tail(i + 2).isZero() || !col2.tail(i + 2).isZero())
                return false;
            ++i;
        }
    }
    return true;
}

template<class MatrixType>
bool realSchurTest(const LValueMatrix<MatrixType>& mat, double precision) {
    RealSchur schur(mat, true);
    if (!isUpperQuasiTriangle(schur.getMatrixT()))
        return false;
    MatrixType A = schur.getMatrixU() * (schur.getMatrixT() * schur.getMatrixU().transpose()).compute();
    if (!matrixNear(A, mat, precision))
        return false;
    return true;
}

int main() {
    {
        using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3>;
        const MatrixType mat1{{-149, 537, -27}, {-50, 180, -9}, {-154, 546, -25}};
        if (!realSchurTest(mat1, 1E-14))
            return 1;
        const MatrixType mat2{{-0.590316, -2.19514, -2.37463},
                              {-1.25006, -0.297493, 1.40349},
                              {0.517063, -0.956614, -0.920775}};
        if (!realSchurTest(mat2, 1E-14))
            return 1;
    }
    return 0;
}
