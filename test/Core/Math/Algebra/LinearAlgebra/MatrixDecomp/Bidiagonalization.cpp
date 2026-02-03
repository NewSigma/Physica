/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Bidiagonalization.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using MatrixType = DenseMatrix<ScalarType, MatrixMajor::Col>;

namespace {
    template<Matrix M>
    void testDecomp(const M& source, double tolerance) noexcept {
        Bidiagonalization obj(source);
        M U = obj.getMatrixU();
        M V = obj.getMatrixV();
        M B = obj.getMatrixB();

        M A = (U * B).compute() * V.transpose();
        expect(matrixNear(A, source, tolerance));
    }
}

int main() {
    {
        const MatrixType mat{{1, 2, 3}, {2, 1, 1}, {-2, 0, 1}};
        testDecomp(mat, 1E-15);
    }
    {
        const MatrixType mat{{1, 2, 3, 4, 5}, {5, 6, 7, 8, 9}, {9, 10, 11, 12, 13}, {7, 6, -8, -9, 5}};
        testDecomp(mat, 1E-14);
    }
    return 0;
}
