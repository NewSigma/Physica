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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/QRDecomp.h"

using namespace Physica;

template<Scalar T, Matrix M>
static void testQR(const QRDecomp<T, false>& qr, const M& m, double prec) noexcept {
    DenseMatrix<T> matrixQ = qr.getMatrixQ();
    M matrixR = qr.getMatrixR();
    M result = matrixQ * matrixR;
    if (!matrixNear(result, m, prec))
        exit(1);
}

template<Scalar T, Matrix M>
static void testQRP(const QRDecomp<T, true>& qr, const M& m, double prec) noexcept {
    DenseMatrix<T> matrixQ = qr.getMatrixQ();
    M matrixR = qr.getMatrixR();
    M result0 = matrixQ * matrixR;
    M result1 = result0 * qr.getMatrixP();
    if (!matrixNear(result1, m, prec))
        exit(1);
}

static void testDecomp(const Matrix auto& m, double prec) noexcept {
    using T = std::remove_cvref<decltype(m)>::type::ScalarType;
    QRDecomp<T, false> qr(m.getRow(), m.getCol());
    QRDecomp<T, true> qrp(m.getRow(), m.getCol());
    if constexpr (HasMKL()) {
        qr.compute_mkl(m);
        qrp.compute_mkl(m);
        testQR(qr, m, prec);
        testQRP(qrp, m, prec);
    }

    qr.compute_base(m);
    testQR(qr, m, prec);
}

int main() {
    {
        using T = float64;
        using Matrix3D = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 3, 3>;
        using Matrix4D = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 4, 4>;
        const Matrix3D m1{2, 3, 4, 1, 1, 9, 1, 2, -6};
        const Matrix4D m2{0, 0.125, 0.125, 0, 0.125, 0, 0, 0.125, 0.125, 0, 0, 0.125, 0, 0.125, 0.125, 0};
        testDecomp(m1, 1E-15);
        testDecomp(m2, 1E-15);
    }
    {
        using T = cfloat32;
        using Matrix43 = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 4, 3>;
        Matrix43 m{
                { 0.49671415,  -0.1382643},
                { 0.64768854,  1.52302986},
                {-0.23415337, -0.23413696},
                { 1.57921282,  0.76743473},
                {-0.46947439, -0.56228743},
                {-0.46341769, -1.91328024},
                {-1.02447039,  1.11792545},
                {-0.17242821, -0.86175486},
                {-0.86175486,  0.34268051},
                { 0.24196227, -1.04525372},
                {-1.11731035,  0.53277921},
                { 1.12050268,  0.50288053}
        };
        testDecomp(m, 1E-6);
    }
    return 0;
}
