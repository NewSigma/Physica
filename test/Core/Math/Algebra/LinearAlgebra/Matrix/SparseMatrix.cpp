/*
 * Copyright 2022-2026 Weibo He.
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
