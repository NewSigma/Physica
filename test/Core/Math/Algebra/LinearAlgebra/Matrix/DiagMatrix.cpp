/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;
using T = float64;
using Tc = cfloat64;

int main() {
    auto matD = DiagMatrix<T>(VectorND<T>::random_uniform<Random<>>(16));
    auto x = VectorND<T>::random_uniform<Random<>>(16);
    VectorND<T> result = hadamard(x, matD.diag());
    x = matD * x;
    expect<RandomSource>(x == result); // Test that there is no aliasing problem
    return 0;
}
