/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void cublas_real() {
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

    void cublas_complex() {
        using T = float32;
        constexpr size_t M = 32;
        constexpr size_t N = 16;
        constexpr size_t K = 24;
        constexpr double Prec = 1E-4;
        using ColMatrix = DenseMatrix<T, MatrixMajor::Col>;
        using RowMatrix = DenseMatrix<T, MatrixMajor::Row>;

        const auto lhs = ColMatrix::template random_uniform<RandomSource>(M, N);
        const auto rhs = ColMatrix::template random_uniform<RandomSource>(N, K);
        const MatrixND<T> answer = lhs * rhs;
        device_obj<ColMatrix> target(M, K);
        (lhs.toDevice() * rhs.toDevice()).assign_cublas(target);
        expect<RandomSource>(matrixNear(target.toHost(), answer, Prec));

        const auto rowLhs = RowMatrix::template random_uniform<RandomSource>(N, M);
        const MatrixND<T> answerT = rowLhs.transpose() * rhs;
        device_obj<ColMatrix> targetT(M, K);
        (rowLhs.toDevice().transpose() * rhs.toDevice()).assign_cublas(targetT);
        expect<RandomSource>(matrixNear(targetT.toHost(), answerT, Prec));

        const auto colLhs = ColMatrix::template random_uniform<RandomSource>(N, M);
        const MatrixND<T> answerC = colLhs.transpose() * rhs;
        device_obj<ColMatrix> targetC(M, K);
        (colLhs.toDevice().transpose() * rhs.toDevice()).assign_cublas(targetC);
        expect<RandomSource>(matrixNear(targetC.toHost(), answerC, Prec));

        const auto initial = ColMatrix::template random_uniform<RandomSource>(M, K);
        device_obj<ColMatrix> accumulated(initial);
        (lhs.toDevice() * rhs.toDevice()).assign_add(accumulated);
        expect<RandomSource>(matrixNear(accumulated.toHost(), MatrixND<T>(initial + answer), Prec));
    }
}

int main() {
    cublas_real();
    cublas_complex();
    return 0;
}
