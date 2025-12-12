/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.cuh"

using namespace Physica;
using T = float32;
using RandomSource = Random<PCG64DXSM, 1000>;

namespace {
    void testLU() {
        auto answer = MatrixND<T>::random_uniform<RandomSource>(32, 32);
        answer.diag() = T(10); // Diagonal dominance for numerical stability
        device_obj<DenseLU<float32, false>> lu(answer.toDevice());

        auto matrixLU = lu.getMatrixLU().toHost();
        MatrixND<T> matrixL = matrixLU.tril();
        matrixL.diag() = float64(1);
        MatrixND<T> result = matrixL * matrixLU.triu();
        if (!matrixNear(result, answer, 1E-3))
            exit(1);

        if (!scalarNear(lu.lnAbsDet(), answer.lnAbsDet(), 1E-3))
            exit(1);

        if (lu.sgndet() != answer.sgndet())
            exit(1);
    }

    void testSolve() {
        const auto A = MatrixND<T>::random_uniform<RandomSource>(32, 32);
        const auto b = MatrixND<T>::random_uniform<RandomSource>(32, 1);
        const VectorND<T> answer = MatrixND<T>(A.inv()) * b.col(0);

        device_obj<DenseLU<float32, true>> lu(A.toDevice());
        auto resultD = b.toDevice();
        lu.solve(resultD);
        const VectorND<T> result = resultD.toHost().col(0);
        if (!vectorNear(result, answer, 1E-3))
            exit(1);
    }
}

int main() {
    testLU();
    testSolve();
    return 0;
}
