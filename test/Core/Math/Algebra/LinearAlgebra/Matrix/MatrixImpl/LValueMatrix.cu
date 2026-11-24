/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using MatrixType = MatrixND<float32>;

namespace {
    void diag() {
        const MatrixType A = MatrixType::random_uniform<Random<>>(4, 3);
        const auto d_A = A.toDeviceAsync();
        device_obj<VectorND<float32>> result = d_A.diag();
        expect(A.diag() == result.toHost());
    }

    void minorDiag() {
        const MatrixType A = MatrixType::random_uniform<Random<>>(4, 4);
        const auto d_A = A.toDeviceAsync();
        int shift = Array<int>{1, 2, 3, -1, -2, -3, 0}.template select<Random<>>();

        device_obj<VectorND<float32>> result = d_A.diag(shift);
        expect(A.diag(shift) == result.toHost());
    }

    void view() {
        auto d_A = MatrixType::random_uniform<Random<>>(5, 6).toDevice();
        static_assert(d_A.row(0).isLValueVector());
        static_assert(d_A.col(0).isLValueVector());
        static_assert(d_A.rows(1, 3).isLValueMatrix());
        static_assert(d_A.topRows(2).isLValueMatrix());
        static_assert(d_A.bottomRows(4).isLValueMatrix());
        static_assert(d_A.cols(1, 3).isLValueMatrix());
        static_assert(d_A.leftCols(2).isLValueMatrix());
        static_assert(d_A.rightCols(4).isLValueMatrix());
        static_assert(d_A.topLeftCorner(2, 3).isLValueMatrix());
        static_assert(d_A.topLeftCorner(2).isLValueMatrix());
        static_assert(d_A.topRightCorner(2, 4).isLValueMatrix());
        static_assert(d_A.bottomLeftCorner(3, 2).isLValueMatrix());
        static_assert(d_A.bottomRightCorner(3, 4).isLValueMatrix());
        static_assert(d_A.bottomRightCorner(2).isLValueMatrix());
        static_assert(d_A.block(1, 3, 2, 3).isLValueMatrix());
    }
}

int main() {
    diag();
    minorDiag();
    view();
    return 0;
}
