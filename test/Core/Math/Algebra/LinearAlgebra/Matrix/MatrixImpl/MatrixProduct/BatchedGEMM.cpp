/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixProduct/BatchedGEMM.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
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

        auto aFirst = a.block(0, row, 0, inner);
        auto bFirst = b.block(0, inner, 0, col);
        auto cFirst = c.block(0, row, 0, col);
        BatchedGEMM<T> batchedGemm(aFirst, bFirst, cFirst, batch);
        for (size_t s = 0; s < batch; ++s)
            batchedGemm.batch(s, {a.block(0, row, s * inner, inner).data_handle(),
                                  b.block(0, inner, s * col, col).data_handle(),
                                  c.block(0, row, s * col, col).data_handle(),});

        if constexpr (HasMKL()) {
            batchedGemm.compute(T(1), T(0));
            expect<RandomSource>(matrixNear(c, answer, 1E-4));

            const std::ptrdiff_t strideA = a.block(0, row, inner, inner).data_handle() - aFirst.data_handle();
            const std::ptrdiff_t strideB = b.block(0, inner, col, col).data_handle() - bFirst.data_handle();
            const std::ptrdiff_t strideC = c.block(0, row, col, col).data_handle() - cFirst.data_handle();
            c.junk();
            batchedGemm.compute(strideA, strideB, strideC, T(1), T(0));
            expect<RandomSource>(matrixNear(c, answer, 1E-4));
        }
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
