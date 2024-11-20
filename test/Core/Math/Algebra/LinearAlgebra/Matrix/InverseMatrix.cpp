/*
 * Copyright 2021-2024 Weibo He.
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
    typedef DenseMatrix<float64, MatrixOption::Col | MatrixOption::Vector, 4, 4> Matrix4x4;
    const Matrix4x4 input{{1, 1, 1, 1}, {1, 1, -1, -1}, {1, -1, 1, -1}, {1, -1, -1, 1}};
    InverseMatrix<Matrix4x4> inv(input);
    const Matrix4x4 result(inv);
    const Matrix4x4 answer{{0.25, 0.25, 0.25, 0.25}, {0.25, 0.25, -0.25, -0.25}, {0.25, -0.25, 0.25, -0.25}, {0.25, -0.25, -0.25, 0.25}};
    if (!matrixNear(answer, result, 1E-15))
        return 1;
    return 0;
}
