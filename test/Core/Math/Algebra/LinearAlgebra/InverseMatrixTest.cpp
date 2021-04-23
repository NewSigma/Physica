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
    typedef DenseMatrix<double, DenseMatrixType::Column | DenseMatrixType::Vector, 3, 3, 3, 3> Matrix4x4;
    Matrix4x4 mat1{{1, 1, 1, 1}, {1, 1, -1, -1}, {1, -1, 1, -1}, {1, -1, -1, 1}};
    InverseMatrix inv(mat1);
    Matrix4x4 result(inv);
    Matrix4x4 answer{{0.25, 0.25, 0.25, 0.25}, {0.25, 0.25, -0.25, -0.25}, {0.25, -0.25, 0.25, -0.25}, {0.25, -0.25, -0.25, 0.25}};
    result -= answer;
    for (const auto& vector : result)
        for (const auto& elemant : vector)
            if (element > 0.001)
                return 1;
    return 0;
}
