/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Test.h"

using namespace Physica;

namespace {
    void testInv(const Matrix auto& m, double prec) {
        using M = std::remove_cvref<decltype(m)>::type;
        using T = M::ScalarType;
        const M result = m.inv();
        M prod = result * m;
        expect(matrixNear(prod, IdentityMatrix<T>(m.getRow()), prec));
    }
}

int main() {
    /* Col major */ {
        using Matrix4D = DenseMatrix<float64, MatrixMajor::Col, 4, 4>;
        testInv(Matrix4D{
                {1,  1,  1,  1},
                {1,  1, -1, -1},
                {1, -1,  1, -1},
                {1, -1, -1,  1}
        }, 1E-14);
        testInv(Matrix4D{
                {1,  2, 0, 1},
                {0,  1, 1, 0},
                {2,  0, 1, 1},
                {1,  1, 0, 1}
        }, 1E-14);
    }
    /* Row major */ {
        using Matrix4D = DenseMatrix<float64, MatrixMajor::Row, 4, 4>;
        testInv(Matrix4D{
                {1,  1,  1,  1},
                {1,  1, -1, -1},
                {1, -1,  1, -1},
                {1, -1, -1,  1}
        }, 1E-14);
        testInv(Matrix4D{
                {1,  2, 0, 1},
                {0,  1, 1, 0},
                {2,  0, 1, 1},
                {1,  1, 0, 1}
        }, 1E-14);
    }
    return 0;
}
