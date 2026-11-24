/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Test.h"

using namespace Physica;

int main() {
    SymmArray<float64> s(3, 0);
    s[1, 0] = 5;
    expect(s[0, 1] == float64(5));

    DenseSymmMatrix<float64> m(3, 0);
    m[1, 0] = 5;
    expect(m[0, 1] == float64(5));
    return 0;
}
