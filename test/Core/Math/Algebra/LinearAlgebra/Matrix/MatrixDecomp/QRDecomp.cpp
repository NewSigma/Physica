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
using Matrix4D = DenseMatrix<float64, MatrixOption::Col | MatrixOption::Element, 4, 4>;

template<Matrix M>
void composite(const QRDecomp<float64>& qr, const M& m) {
    M matrixQ = qr.getMatrixQ();
    M matrixR = qr.getMatrixR();
    M result = matrixQ * matrixR;
    if (!matrixNear(result, m, 1E-15))
        exit(1);
}

template<Matrix M>
void testQR(const M& m) {
    QRDecomp<float64> qr(m.getRow(), m.getCol());
    if constexpr (HasMKL()) {
        qr.compute_mkl(m);
        composite(qr, m);
    }
    qr.compute_base(m);
    composite(qr, m);
}

int main() {
    const Matrix3D m1{2, 3, 4, 1, 1, 9, 1, 2, -6};
    const Matrix4D m2{0, 0.125, 0.125, 0, 0.125, 0, 0, 0.125, 0.125, 0, 0, 0.125, 0, 0.125, 0.125, 0};
    testQR(m1);
    testQR(m2);
    return 0;
}
