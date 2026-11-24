/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Cholesky.h"
#include "Test.h"

using namespace Physica;

int main() {
    {
        using Matrix3x3 = DenseMatrix<float64, MatrixMajor::Col, 3, 3>;
        Matrix3x3 mat{5, -2, 0, -2, 3, -1, 0, -1, 1};
        Cholesky cholesky(mat);
        Matrix3x3 decomp(cholesky);
        Matrix3x3 answer{sqrt(5), -2 / sqrt(5), 0, 0, sqrt(11.0 / 5), -sqrt(5.0 / 11), 0, 0, sqrt(6.0 / 11)};
        decomp -= answer;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                expect(decomp[i, j] < 0.001);
    }
    {
        using Matrix3x3 = DenseMatrix<float64, MatrixMajor::Col, 3, 3>;
        Matrix3x3 mat{{5, -2, 0}, {-2, 3, -1}, {0, -1, 1}};
        Cholesky cholesky(mat);
        Matrix3x3 decomp(cholesky);
        Matrix3x3 answer{{sqrt(5), -2 / sqrt(5), 0}, {0, sqrt(11.0 / 5), -sqrt(5.0 / 11)}, {0, 0, sqrt(6.0 / 11)}};
        decomp -= answer;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                expect(decomp[i, j] < 0.001);
    }
    {
        using Matrix5x5 = DenseMatrix<float64, MatrixMajor::Col, 5, 5>;
        const Matrix5x5 mat{
            {1.066666667,             0,  0.1523809524,             0, 0.05079365079},
            {            0,  0.1523809524,             0, 0.05079365079,             0},
            {0.1523809524,             0, 0.05079365079,             0, 0.02308802309},
            {            0, 0.05079365079,             0, 0.02308802309,             0},
            {0.05079365079,             0, 0.02308802309,             0, 0.01243201243}
        };
        const Matrix5x5 answer{
            {{1.0328, 0., 0., 0., 0.},
            {0., 0.39036, 0., 0., 0.},
            {0.147542, 0., 0.170367, 0., 0.},
            {0., 0.13012, 0., 0.0784653, 0.},
            {0.0491807, 0., 0.0929275, 0., 0.037118}}
        };
        const Matrix5x5 result = Cholesky(mat);
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                expect(scalarNear(result[i, j], answer[j, i], 1E-5));
    }
    return 0;
}
