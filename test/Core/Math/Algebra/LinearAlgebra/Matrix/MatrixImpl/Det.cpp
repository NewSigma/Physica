/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica;
using Matrix4D = DenseMatrix<float64, MatrixMajor::Col, 4, 4>;

int main() {
    const Matrix4D m2{0, 0.125, 0.125, 0, 0.125, 0, 0, 0.125, 0.125, 0, 0, 0.125, 0, 0.125, 0.125, 0};
    Matrix4D e2 = (exp(m2));
    return !scalarNear(float64(1), e2.det(), 1E-15);
}
