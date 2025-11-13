/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"

using namespace Physica;
using T = float64;
using Matrix3D = DenseMatrix<T, MatrixOption::Col, 3, 3>;

int main() {
    const Matrix3D answer{1, 100,  10000, 0.01, 1, 100, 0.0001, 0.01, 1};
    Matrix3D m = answer;
    const auto D = DiagMatrix<T>(m.balance());
    Matrix3D result = D.inv() * m * D;
    if (!matrixNear(result, answer, 1E-15))
        return 1;
    return 0;
}
