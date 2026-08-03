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
