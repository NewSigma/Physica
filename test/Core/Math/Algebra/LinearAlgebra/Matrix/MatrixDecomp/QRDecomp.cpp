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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/QRDecomp.h"

using namespace Physica;
using Matrix3D = DenseMatrix<float64, MatrixOption::Col | MatrixOption::Element, 3, 3>;

void testQR(const QRDecomp<float64>& qr, const Matrix3D& answer) {
    Matrix3D matrixQ = qr.getMatrixQ();
    Matrix3D matrixR = qr.getMatrixR();
    Matrix3D result = matrixQ * matrixR;
    if (!matrixNear(result, answer, 1E-15))
        exit(1);
}

int main() {
    const Matrix3D answer{2, 3, 4, 1, 1, 9, 1, 2, -6};
    QRDecomp<float64> qr(3, 3);
    if constexpr (HasMKL()) {
        qr.compute_mkl(answer);
        testQR(qr, answer);
    }
    qr.compute_base(answer);
    testQR(qr, answer);
    return 0;
}
