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
    {
        typedef DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Element, 3, 3, 3, 3> Matrix3x3;
        Matrix3x3 mat{5, -2, 0, -2, 3, -1, 0, -1, 1};
        Cholesky cholesky(mat);
        Matrix3x3 decomp(cholesky);
        Matrix3x3 answer{sqrt(5), -2 / sqrt(5), 0, 0, sqrt(11.0 / 5), -sqrt(5.0 / 11), 0, 0, sqrt(6.0 / 11)};
        decomp -= answer;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (decomp(i, j) > 0.001)
                    return 1;
    }
    {
        typedef DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3, 3, 3> Matrix3x3;
        Matrix3x3 mat{{5, -2, 0}, {-2, 3, -1}, {0, -1, 1}};
        Cholesky cholesky(mat);
        Matrix3x3 decomp(cholesky);
        Matrix3x3 answer{{sqrt(5), -2 / sqrt(5), 0}, {0, sqrt(11.0 / 5), -sqrt(5.0 / 11)}, {0, 0, sqrt(6.0 / 11)}};
        decomp -= answer;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (decomp(i, j) > 0.001)
                    return 1;
    }
    return 0;
}
