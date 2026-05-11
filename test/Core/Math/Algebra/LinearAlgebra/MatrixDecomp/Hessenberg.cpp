/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Hessenberg.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Test.h"

using namespace Physica;

namespace {
    bool isHessenbergMatrix(const Matrix auto& m) {
        const size_t order = m.getRow();
        if (m.getRow() != m.getCol())
            return false;
        if (order <= 2)
            return true;
        for (size_t i = 0; i < order - 2; ++i) {
            for (size_t j = i + 2; j < order; ++j)
                if (!m[j, i].isZero())
                    return false;
        }
        return true;
    }

    template<Matrix M>
    void hessTest(const M& source, double tolerance) {
        Hessenberg<typename M::ScalarType> hess(source);
        M H = hess.getMatrixH();
        expect(isHessenbergMatrix(H));

        M Q = hess.getMatrixQ();
        M A = (Q * H).compute() * Q.hermite();
        expect(matrixNear(A, source, tolerance));
    }
}

int main() {
    using RealType = float64;
    using ComplexType = Complex<RealType>;
    {
        using MatrixType = DenseMatrix<RealType, MatrixMajor::Col, 4, 4>;
        const MatrixType mat{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        hessTest(mat, 1E-15);
    }
    {
        using MatrixType = DenseMatrix<RealType, MatrixMajor::Col, 3, 3>;
        const MatrixType mat{{-149, 537, -27}, {-50, 180, -9}, {-154, 546, -25}};
        hessTest(mat, 1E-15);
    }
    /* Test degeneracy */ {
        using MatrixType = DenseMatrix<RealType, MatrixMajor::Col, 6, 6>;
        const MatrixType mat1{{ 0.1343184046,             0,             0, -0.1343184056,             0,             0},
                              {            0,  0.1341424528,             0,             0, -0.1341424541,             0},
                              {            0,             0,  0.1342191829,             0,             0, -0.1342191848},
                              {-0.1343184056,             0,             0,  0.1343184065,             0,             0},
                              {            0, -0.1341424541,             0,             0,  0.1341424554,             0},
                              {            0,             0, -0.1342191848,             0,             0,  0.1342191868}};
        hessTest(mat1, 1E-15);

        using ComplexMatrix = DenseMatrix<ComplexType, MatrixMajor::Col, 6, 6>;
        const ComplexMatrix mat2 = mat1;
        hessTest(mat2, 1E-15);
    }
    {
        using MatrixType = DenseMatrix<RealType, MatrixMajor::Col, 8, 8>;
        const MatrixType mat{{2.5, 0, 0.866025, 0, 0.481812, 0, 0.318105, 0},
                            {0, 10.5, 0, 4.97494, 0, 3.06186, 0, 2.12934},
                            {0.866025, 0, 25.5, 0, 14.1869, 0, 9.36654, 0},
                            {0, 4.97494, 0, 49.5, 0, 30.4651, 0, 21.1867},
                            {0.481812, 0, 14.1869, 0, 84.5, 0, 55.789, 0},
                            {0, 3.06186, 0, 30.4651, 0, 132.5, 0, 92.1458},
                            {0.318105, 0, 9.36654, 0, 55.789, 0, 195.5, 0},
                            {0, 2.12934, 0, 21.1867, 0, 92.1458, 0, 275.5}};
        hessTest(mat, 1E-13);
    }
    /* Complex case */ {
        using MatrixType = DenseMatrix<ComplexType, MatrixMajor::Col, 3, 3>;
        const MatrixType mat{{{2, 1}, {-3, 6}, {12, 7}}, {{-50, -9}, {2, 180}, {-9, -6}}, {{-7, 8}, {546, 0}, {0, -25}}};
        hessTest(mat, 1E-12);
    }
    return 0;
}
