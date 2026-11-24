/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/SparseMatrix.h"
#include "Test.h"

using namespace Physica;
using T = float64;

static_assert(SparseMatrix<T>::isSparse());
static_assert(!DenseMatrix<T>::isSparse());

int main() {
    {
        SparseMatrix<T> mat(3, 3);
        mat.insert(1, 0, 1);
        mat.insert(1, 1, 0);
        mat.insert(1, 1, 2);
        mat.insert(1, 2, 1);
        {
            DenseMatrix<T> answer{
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            };
            expect(matrixNear(mat, answer, 1E-16));
        }
        {
            VectorND<T> v{1, 1, 1};
            VectorND<T> result = mat * v;
            VectorND<T> answer{1, 2, 1};
            expect(vectorNear(result, answer, 1E-16));
        }
    }
    {
        SparseMatrix<T> mat(4, 4);
        mat.insert(1, 0, 0);
        mat.insert(4, 1, 1);
        mat.insert(2, 1, 0);
        mat.insert(3, 0, 1);
        DenseMatrix<T> answer{
            {1, 2, 0, 0},
            {3, 4, 0, 0},
            {0, 0, 0, 0},
            {0, 0, 0, 0}
        };
        expect(matrixNear(mat, answer, 1E-15));
    }
    return 0;
}
