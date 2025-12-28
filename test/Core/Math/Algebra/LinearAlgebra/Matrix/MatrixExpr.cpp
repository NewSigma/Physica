/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica;

int main() {
    {
        DenseMatrix<float64, MatrixOption::Col, 3, 3> mat1{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
        DenseMatrix<float32, MatrixOption::Col, 3, 3> mat2{1, 1, 1, 1, 1, 1, 1, 1, 1};
        {
            DenseMatrix<float64, MatrixOption::Row, 3, 3> mat = -(mat1 + mat2);
            for (size_t i = 0; i < mat.getRow(); ++i)
                for (size_t j = 0; j < mat.getCol(); ++j)
                    if (mat[i, j] != float64(-2))
                        return 1;
        }
        {
            DenseMatrix<float64, MatrixOption::Row, 3, 3> mat = mat1 * mat2;
            for (size_t i = 0; i < mat.getRow(); ++i)
                for (size_t j = 0; j < mat.getCol(); ++j)
                    if (mat[i, j] != float64(3))
                        return 1;
        }
    }
    /* ContinuousMatrixBlock<Derived, 1, 1> */ {
        using ScalarType = float64;
        DenseMatrix<ScalarType, MatrixOption::Col, Physica::Dynamic, 1> mat(2, 1);
        mat[0, 0] = 1.0;
        mat[1, 0] = 2.0;
        if (mat.row(1)[0] != ScalarType(2))
            return 1;
    }
    return 0;
}
