/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void assign() {
        constexpr int N = 16;
        using T = float64;
        const auto a = VectorND<T>::random_uniform<RandomSource>(N);
        const auto b = VectorND<T>::random_uniform<RandomSource>(N);
        auto y = MatrixND<T>::junk(N, N);
        y = a * b.transpose();
        expect<RandomSource>(matrixNear(y, MatrixND<T>(a * b.transpose()), 1E-15));
    }

    void sum() {
        constexpr int N = 16;
        using T = float64;
        const auto a = VectorND<T>::random_uniform<RandomSource>(N);
        const auto b = VectorND<T>::random_uniform<RandomSource>(N);
        const MatrixND<T> y = a * b.transpose();
        expect<RandomSource>(scalarNear((a * b.transpose()).sum(), y.sum(), 1E-15));
    }
}

int main() {
    sum();
    assign();
    return 0;
}
