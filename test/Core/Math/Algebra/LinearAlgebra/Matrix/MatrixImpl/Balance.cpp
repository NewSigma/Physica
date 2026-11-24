/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Matrix3D = DenseMatrix<T, MatrixMajor::Col, 3, 3>;

int main() {
    const Matrix3D answer{1, 100,  10000, 0.01, 1, 100, 0.0001, 0.01, 1};
    Matrix3D m = answer;
    const auto D = DiagMatrix<T>(m.balance());
    Matrix3D result = D.inv() * m * D;
    expect(matrixNear(result, answer, 1E-15));
    return 0;
}
