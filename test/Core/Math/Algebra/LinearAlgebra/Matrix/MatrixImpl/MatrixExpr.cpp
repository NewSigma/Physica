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
#include "Test.h"

using namespace Physica;

namespace {
    void general() {
        DenseMatrix<float64, MatrixMajor::Col, 3, 3> mat1{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
        DenseMatrix<float32, MatrixMajor::Col, 3, 3> mat2{1, 1, 1, 1, 1, 1, 1, 1, 1};
        {
            DenseMatrix<float64, MatrixMajor::Row, 3, 3> mat = -(mat1 + mat2);
            for (size_t i = 0; i < mat.getRow(); ++i)
                for (size_t j = 0; j < mat.getCol(); ++j)
                    expect(mat[i, j] == float64(-2));
        }
        {
            DenseMatrix<float64, MatrixMajor::Row, 3, 3> mat = mat1 * mat2;
            for (size_t i = 0; i < mat.getRow(); ++i)
                for (size_t j = 0; j < mat.getCol(); ++j)
                    expect(mat[i, j] == float64(3));
        }
    }

    void block1x1() {
        using ScalarType = float64;
        DenseMatrix<ScalarType, MatrixMajor::Col, Physica::Dynamic, 1> mat(2, 1);
        mat[0, 0] = 1.0;
        mat[1, 0] = 2.0;
        expect(mat.row(1)[0] == ScalarType(2));
    }
    /**
     * Test that we are free of a set of lifetime problems under CXX23
     */
    void lifetimeCXX23() {
        using T = float32;
        auto diag = (MatrixND<T>::identity(3) + MatrixND<T>(3, 3, 5)).diag();
        for (int i = 0; i < 3; ++i)
            expect(diag.calc(i) == T(6));
    }
}

int main() {
    general();
    block1x1();
    lifetimeCXX23();
    return 0;
}
