/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseQR.cuh"

using namespace Physica;
using Matrix3D = DenseMatrix<float32, MatrixMajor::Col, 3>;

namespace {
    void testQR(device_obj<DenseQR<float32>>& qr, const Matrix3D& answer) {
        device_obj<Matrix3D> matrixQ = qr.getMatrixQ();
        device_obj<Matrix3D> matrixR = qr.getMatrixR();
        Matrix3D result = device_obj<Matrix3D>(matrixQ * matrixR).toHost();
        if (!matrixNear(result, answer, 1E-6))
            exit(1);
    }
}

int main() {
    const Matrix3D answer{2, 3, 4, 1, 1, 9, 1, 2, -6};
    device_obj<DenseQR<float32>> qr(3, 3);

    qr.compute(answer.toDevice());
    testQR(qr, answer);
    return 0;
}
