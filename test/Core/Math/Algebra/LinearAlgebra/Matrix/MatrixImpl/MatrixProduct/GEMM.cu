/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void cublas() {
        using T = float64;
        constexpr size_t M = 64;
        constexpr size_t N = 32;
        constexpr size_t K = 16;
        using ColMatrix = DenseMatrix<T, MatrixMajor::Col>;
        using RowMatrix = DenseMatrix<T, MatrixMajor::Row>;

        const auto colMatrix = ColMatrix::random_uniform<RandomSource>(M, N);
        {
            const auto rhs = ColMatrix::random_uniform<RandomSource>(N, K);
            const MatrixND<T> answer = colMatrix * rhs;
            device_obj<ColMatrix> target(M, K);
            (colMatrix.toDevice() * rhs.toDevice()).assign_cublas(target);
            expect<RandomSource>(matrixNear(target.toHost(), answer, 1E-13));
        }
        {
            const auto rowMatrix = RowMatrix::random_uniform<RandomSource>(M, N);
            const auto rhs = ColMatrix::random_uniform<RandomSource>(M, K);
            const MatrixND<T> answer = rowMatrix.transpose() * rhs;
            device_obj<ColMatrix> target(N, K);
            (rowMatrix.toDevice().transpose() * rhs.toDevice()).assign_cublas(target);
            expect<RandomSource>(matrixNear(target.toHost(), answer, 1E-13));
        }
    }
}

int main() {
    cublas();
    return 0;
}
