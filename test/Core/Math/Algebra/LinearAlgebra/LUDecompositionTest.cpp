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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;

int main() {
    typedef DenseMatrix<Scalar<Double, false>, DenseMatrixType::Row | DenseMatrixType::Vector, 3, 3, 3, 3> Matrix3x3;
    Matrix3x3 mat1{{2, 3, 4}, {1, 1, 9}, {1, 2, -6}};
    LUDecomposition lu(mat1);
    Matrix3x3 decomp(lu);
    Matrix3x3 answer{{2, 3, 4}, {0.5, -0.5, 7}, {0.5, -1, -1}};
    decomp -= answer;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (decomp(i, j) > 0.001)
                return 1;
    return 0;
}
