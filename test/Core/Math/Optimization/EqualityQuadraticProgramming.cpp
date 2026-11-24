/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2022 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Optimization/QuadraticProgramming/EqualityQuadraticProgramming.h"
#include "Test.h"

using namespace Physica;

int main() {
    using ScalarType = float64;
    {
        using VectorType = Vector3D<ScalarType>;
        const DenseMatrix<ScalarType, MatrixMajor::Row, 3, 3> G{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const VectorType c{0, 0, 0};
        const DenseMatrix<ScalarType> constraints{};
        const VectorType answer{0, 0, 0};

        const VectorType initial{1, -1, 1};
        EqualityQuadraticProgramming<ScalarType> QP(G, c, constraints, initial);
        expect(vectorNear(QP.getSolution(), answer, 1E-16));
    }
    {
        using VectorType = Vector3D<ScalarType>;
        const DenseMatrix<ScalarType, MatrixMajor::Row, 3, 3> G{{6, 2, 1}, {2, 5, 2}, {1, 2, 4}};
        const VectorType c{-8, -3, -3};
        const DenseMatrix<ScalarType, MatrixMajor::Row, 2, 4> constraints{{1, 0, 1, 3}, {0, 1, 1, 0}};
        const VectorType answer{2, -1, 1};

        const VectorType initial{1, 1, 1};
        EqualityQuadraticProgramming<ScalarType> QP(G, c, constraints, initial);
        expect(vectorNear(QP.getSolution(), answer, 1E-15)); //Precision can reach 1E-16 in debug mode
    }
    return 0;
}
