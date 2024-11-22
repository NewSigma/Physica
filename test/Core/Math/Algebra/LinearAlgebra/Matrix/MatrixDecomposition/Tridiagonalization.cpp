/*
 * Copyright 2022-2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Tridiagonalization.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica::Core;

template<Matrix M>
bool doTest(const M& source, double tolerance) {
    using ScalarType = M::ScalarType;
    using GeneralMatrix = DenseMatrix<ScalarType>;
    Tridiagonalization<ScalarType> tri(source);
    GeneralMatrix T = tri.getMatrixT();
    GeneralMatrix Q = tri.getMatrixQ();
    GeneralMatrix A = (Q * T).compute() * Q.hermite();
    if (!matrixNear(A, source, tolerance))
        return false;
    return true;
}

int main() {
    using RealType = float64;
    {
        using MatrixType = DenseMatrix<RealType, MatrixOption::Col | MatrixOption::Vector, 4, 4>;
        const MatrixType temp{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        const MatrixType mat = temp + temp.transpose();
        if (!doTest(mat, 1E-14))
            return 1;
    }
    {
        using MatrixType = DenseSymmMatrix<RealType>;
        const auto mat = MatrixType::random_uniform(8, Random<MT19937, std::mt19937::default_seed>::getInstance());
        if (!doTest(mat, 1E-12))
            return 1;
    }
    /* Complex case */ {
        using MatrixType = DenseMatrix<Complex<RealType>, MatrixOption::Col | MatrixOption::Vector, 3, 3>;
        const MatrixType temp{{{2, 1}, {-3, 6}, {12, 7}}, {{-50, -9}, {2, 180}, {-9, -6}}, {{-7, 8}, {546, 0}, {0, -25}}};
        const MatrixType mat = temp + temp.hermite();
        if (!doTest(mat, 1E-12))
            return 1;
    }
    return 0;
}
