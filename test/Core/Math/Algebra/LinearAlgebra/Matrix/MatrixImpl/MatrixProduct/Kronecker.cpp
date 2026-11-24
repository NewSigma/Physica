/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MCG>;
using T = float64;

int main() {
    const auto a = MatrixND<T>::random_uniform<RandomSource>(4, 4);
    const auto b = MatrixND<T>::random_uniform<RandomSource>(6, 6);
    auto k = MatrixND<T>::junk(a.getRow() * b.getRow(), a.getCol() * b.getCol());
    k = kronecker(a, b);
    expect<RandomSource>(matrixNear(k, MatrixND<T>(kronecker(a, b)), 1E-15));

    const auto v = VectorND<T>::random_uniform<RandomSource>(kronecker(a, b).getRow());
    expect<RandomSource>(vectorNear(k * v, kronecker(a, b) * v, 4UL));
    return 0;
}
