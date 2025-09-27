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

namespace {
    template<Scalar T, Matrix M, bool Pivot>
    void testQR(const QRDecomp<T, Pivot>& qr, const M& m, double decompPrec, double detPrec) noexcept {
        DenseMatrix<T> matrixQ = qr.getMatrixQ();
        M matrixR = qr.getMatrixR();
        M result = matrixQ * matrixR;
        if constexpr (Pivot)
            result = M(result * qr.getMatrixP());

        if (!matrixNear(result, m, decompPrec))
            exit(1);

        if (m.isSquare() && (m.getRow() <= 3) && !scalarNear(m.det(), qr.det(), detPrec))
            exit(1);
    }

    void testDecomp(const Matrix auto& m, double decompPrec, double detPrec) noexcept {
        using T = std::remove_cvref<decltype(m)>::type::ScalarType;
        QRDecomp<T, false> qr(m.getRow(), m.getCol());
        QRDecomp<T, true> qrp(m.getRow(), m.getCol());
        if constexpr (HasMKL()) {
            qr.compute_mkl(m);
            qrp.compute_mkl(m);
            testQR(qr, m, decompPrec, detPrec);
            testQR(qrp, m, decompPrec, detPrec);
        }

        qr.compute_base(m);
        testQR(qr, m, decompPrec, detPrec);
    }
}

int main() {
    {
        using T = float64;
        using Matrix3D = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 3, 3>;
        using Matrix4D = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 4, 4>;
        const Matrix3D m1{2, 3, 4, 1, 1, 9, 1, 2, -6};
        const Matrix4D m2{0, 0.125, 0.125, 0, 0.125, 0, 0, 0.125, 0.125, 0, 0, 0.125, 0, 0.125, 0.125, 0};
        testDecomp(m1, 1E-15, 1E-13);
        testDecomp(m2, 1E-15, 1E-13);
    }
    {
        using T = cfloat64;
        using Matrix3D = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, 3, 3>;
        const Matrix3D m1{
                {0.314168, 0.121569},
                {0.542236, 0.789234},
                {0.681570, 0.478108},
                {0.912647, 0.227165},
                {0.216599, 0.948223},
                {0.006347, 0.337121},
                {0.632359, 0.532767},
                {0.253040, 0.873016},
                {0.218257, 0.103735}
        };
        testDecomp(m1, 1E-14, 1E-14);
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
        testDecomp(m, 1E-6, 0);
    }
    return 0;
}
