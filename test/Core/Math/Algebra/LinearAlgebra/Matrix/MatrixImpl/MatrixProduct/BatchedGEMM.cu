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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixProduct/BatchedGEMM.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    template<int MajorA, int MajorB, int MajorC>
    void gemm() {
        using T = float32;
        const size_t batch = 3, row = 7, col = 5, inner = 6;
        using MatA = DenseMatrix<T, MajorA>;
        using MatB = DenseMatrix<T, MajorB>;
        using MatC = DenseMatrix<T, MajorC>;
        auto a = MatA::template random_uniform<RandomSource>(row, batch * inner);
        auto b = MatB::template random_uniform<RandomSource>(inner, batch * col);
        auto c = MatC::template random_uniform<RandomSource>(row, batch * col);

        MatC answer(row, batch * col);
        for (size_t s = 0; s < batch; ++s) {
            auto answerBlock = answer.block(0, row, s * col, col);
            answerBlock = a.block(0, row, s * inner, inner) * b.block(0, inner, s * col, col);
        }

        device_obj<MatA> aDevice = a.toDevice();
        device_obj<MatB> bDevice = b.toDevice();
        device_obj<MatC> cDevice = c.toDevice();
        auto aFirst = aDevice.block(0, row, 0, inner);
        auto bFirst = bDevice.block(0, inner, 0, col);
        auto cFirst = cDevice.block(0, row, 0, col);
        device_obj<BatchedGEMM<T>> batchedGemm(aFirst, bFirst, cFirst, batch);
        for (size_t s = 0; s < batch; ++s)
            batchedGemm.batch(s, {aDevice.block(0, row, s * inner, inner).data_handle(),
                                  bDevice.block(0, inner, s * col, col).data_handle(),
                                  cDevice.block(0, row, s * col, col).data_handle(),});
        batchedGemm.compute(T(1), T(0));
        expect<RandomSource>(matrixNear(cDevice.toHost(), answer, 1E-4));
    }
}

int main() {
    gemm<MatrixMajor::Row, MatrixMajor::Row, MatrixMajor::Row>();
    gemm<MatrixMajor::Row, MatrixMajor::Row, MatrixMajor::Col>();
    gemm<MatrixMajor::Row, MatrixMajor::Col, MatrixMajor::Row>();
    gemm<MatrixMajor::Row, MatrixMajor::Col, MatrixMajor::Col>();
    gemm<MatrixMajor::Col, MatrixMajor::Row, MatrixMajor::Row>();
    gemm<MatrixMajor::Col, MatrixMajor::Row, MatrixMajor::Col>();
    gemm<MatrixMajor::Col, MatrixMajor::Col, MatrixMajor::Row>();
    gemm<MatrixMajor::Col, MatrixMajor::Col, MatrixMajor::Col>();
    return 0;
}
