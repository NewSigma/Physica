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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/SVD.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    template<Matrix M>
    void decomp(const M& source, double tolerance) {
        SVD<T> svd(source);
        const auto& U = svd.getMatrixU();
        const auto& V = svd.getMatrixV();
        const auto& v = svd.getSignSingulars();

        M A(source.getRow(), source.getCol(), 0);
        for (size_t i = 0; i < v.getLength(); ++i)
            A += U.col(i) * V.col(i).transpose() * v[i];
        expect(matrixNear(A, source, tolerance));
    }
}

int main() {
    {
        MatrixND<T> mat{{1, 2, 3}, {2, 1, 1}, {-2, 0, 1}};
        decomp(mat, 1E-14);
    }
    {
        const MatrixND<T> mat{{1, 2, 3, 4, 5}, {5, 6, 7, 8, 9}, {9, 10, 11, 12, 13}, {7, 6, -8, -9, 5}};
        decomp(mat, 1E-14);
    }
    {
        const MatrixND<T> mat{{1, 0, 0}, {0, -1, 0}, {0, 0, 1}};
        decomp(mat, 1E-14);
    }
    return 0;
}
