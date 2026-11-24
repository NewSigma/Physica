/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using MatrixType = MatrixND<float32>;

namespace {
    void diag() {
        const MatrixType A = MatrixType::random_uniform<Random<>>(4, 3);
        const auto d_A = A.toDeviceAsync();
        device_obj<VectorND<float32>> result = (-d_A).diag();
        expect((-A).diag() == result.toHost());
    }

    void minorDiag() {
        const MatrixType A = MatrixType::random_uniform<Random<>>(4, 4);
        const auto d_A = A.toDeviceAsync();
        int shift = Array<int>{1, 2, 3, -1, -2, -3, 0}.template select<Random<>>();

        device_obj<VectorND<float32>> result = (-d_A).diag(shift);
        expect((-A).diag(shift) == result.toHost());
    }

    void sum() {
        const MatrixType A = MatrixType::random_uniform<Random<>>(16, 16);
        const auto d_A = A.toDeviceAsync();
        expect(scalarNear(d_A.diag().sum(), A.diag().sum(), 1E-6));
    }
}

static_assert(std::same_as<decltype(device_obj<MatrixND<float32>>{}.transpose().transpose()), device_obj<MatrixND<float32>>&&>, "[Error]: transpose-transpose should yield it self");

int main() {
    diag();
    minorDiag();
    sum();
    return 0;
}
