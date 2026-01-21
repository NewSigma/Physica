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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica;
using Matrix4D = DenseMatrix<float64, MatrixOption::Col, 4, 4>;

int main() {
    const Matrix4D m2{0, 0.125, 0.125, 0, 0.125, 0, 0, 0.125, 0.125, 0, 0, 0.125, 0, 0.125, 0.125, 0};
    Matrix4D e2 = (exp(m2));
    return !scalarNear(float64(1), e2.det(), 1E-15);
}
