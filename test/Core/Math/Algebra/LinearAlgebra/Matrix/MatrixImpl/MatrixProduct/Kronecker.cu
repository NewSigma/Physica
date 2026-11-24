/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MCG>;
using T = float32;

int main() {
    const auto a = MatrixND<T>::random_uniform<RandomSource>(2, 2);
    const auto b = MatrixND<T>::random_uniform<RandomSource>(32, 32);
    MatrixND<T> answer = kronecker(a, b);

    auto d_a = a.toDeviceAsync();
    auto d_b = b.toDeviceAsync();
    device_obj<MatrixND<T>> result = kronecker(d_a, d_b);
    expect(answer == result.toHost());

    result += kronecker(d_a, d_b);
    answer *= T(2);
    expect(matrixNear(answer, result.toHost(), 1E-5));
    return 0;
}
