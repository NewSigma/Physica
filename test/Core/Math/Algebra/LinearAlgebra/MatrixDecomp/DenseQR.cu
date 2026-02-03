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
