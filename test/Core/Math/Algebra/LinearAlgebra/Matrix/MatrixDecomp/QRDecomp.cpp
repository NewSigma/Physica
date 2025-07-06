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
using T = float64;
using Matrix3D = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 3, 3>;
using Matrix4D = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 4, 4>;

template<Matrix M>
void testQR(const QRDecomp<T, false>& qr, const M& m) {
    M matrixQ = qr.getMatrixQ();
    M matrixR = qr.getMatrixR();
    M result = matrixQ * matrixR;
    if (!matrixNear(result, m, 1E-15))
        exit(1);
}

template<Matrix M>
void testQRP(const QRDecomp<T, true>& qr, const M& m) {
    M matrixQ = qr.getMatrixQ();
    M matrixR = qr.getMatrixR();
    M result0 = matrixQ * matrixR;
    M result1 = result0 * qr.getMatrixP();
    if (!matrixNear(result1, m, 1E-15))
        exit(1);
}

void testDecomp(const Matrix auto& m) {
    QRDecomp<T, false> qr(m.getRow(), m.getCol());
    QRDecomp<T, true> qrp(m.getRow(), m.getCol());
    if constexpr (HasMKL()) {
        qr.compute_mkl(m);
        qrp.compute_mkl(m);
        testQR(qr, m);
        testQRP(qrp, m);
    }

    qr.compute_mkl(m);
    testQR(qr, m);
}

int main() {
    const Matrix3D m1{2, 3, 4, 1, 1, 9, 1, 2, -6};
    const Matrix4D m2{0, 0.125, 0.125, 0, 0.125, 0, 0, 0.125, 0.125, 0, 0, 0.125, 0, 0.125, 0.125, 0};
    testDecomp(m1);
    testDecomp(m2);
    return 0;
}
