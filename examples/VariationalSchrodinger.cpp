/*
 * Copyright 2021 WeiBo He.
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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Cholesky.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/RealSchur.h"

using namespace Physica::Core;

constexpr size_t baseSetCount = 5;

using ScalarType = Scalar<Double, false>;
using MatrixType = DenseMatrix<ScalarType,
                               DenseMatrixOption::Column | DenseMatrixOption::Vector,
                               baseSetCount,
                               baseSetCount>;

MatrixType getHamiltonMatrix() {
    MatrixType result = MatrixType::Zeros(baseSetCount, baseSetCount);
    for (size_t i = 0; i < baseSetCount; ++i) {
        for (size_t j = i % 2; j < baseSetCount; j += 2) {
            const ScalarType sum = ScalarType(i + j);
            const ScalarType pro = ScalarType(i * j);
            const ScalarType numerator = ScalarType::One() - sum - ScalarType::Two() * pro;
            const ScalarType denominator = (sum + ScalarType(3)) * (sum + ScalarType::One()) * (sum - ScalarType::One());
            result(i, j) = ScalarType(-8) * numerator * denominator;
        }
    }
    return result;
}

MatrixType getOverlapMatrix() {
    MatrixType result = MatrixType::Zeros(baseSetCount, baseSetCount);
    for (size_t i = 0; i < baseSetCount; ++i) {
        for (size_t j = i % 2; j < baseSetCount; j += 2) {
            const ScalarType sum = ScalarType(i + j);
            const ScalarType term1 = ScalarType::Two() * (reciprocal(sum + ScalarType::One()) + reciprocal(sum + ScalarType(5)));
            const ScalarType term2 = ScalarType(4) * reciprocal(sum + ScalarType(3));
            result(i, j) = term1 - term2;
        }
    }
    return result;
}
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:29-34
 */
int main() {
    MatrixType overlap = getOverlapMatrix();
    //Workaround for generalised eigenvalue problem
    MatrixType cholesky = Cholesky(overlap);
    MatrixType inv_cholesky = cholesky.inverse();
    MatrixType hamilton = getHamiltonMatrix();
    MatrixType hamilton_mod = (inv_cholesky.transpose() * hamilton).compute() * inv_cholesky;
    MatrixType schur = RealSchur(hamilton_mod);
    for (size_t i = 0; i < baseSetCount; ++i)
        std::cout << schur(i, i) << '\n';
    return 0;
}
