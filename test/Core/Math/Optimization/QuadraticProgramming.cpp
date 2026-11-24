/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2022 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Optimization/QuadraticProgramming/QuadraticProgramming.h"
#include "Test.h"

using namespace Physica;

int main() {
    using ScalarType = float64;
    {
        using VectorType = Vector2D<ScalarType>;
        const DenseMatrix<ScalarType, MatrixMajor::Row, 2, 2> G{{2, 0}, {0, 2}};
        const VectorType c{-2, -5};
        const DenseMatrix<ScalarType, MatrixMajor::Row> equality{};
        const DenseMatrix<ScalarType, MatrixMajor::Row, 5, 3> inequality{{1, -2, -2}, {-1, -2, -6}, {-1, 2, -2}, {1, 0, 0}, {0, 1, 0}};
        const VectorType answer{1.4, 1.7};

        const VectorType initial{2, 0};
        QuadraticProgramming<ScalarType> QP(G, c, equality, inequality, initial);
        QP.compute();
        expect(vectorNear(QP.getSolution(), answer, 1E-15));
    }
    return 0;
}
