/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;

int main() {
    const VectorND<ScalarType> a{1, 2, 3, 4, 5};
    expect(scalarNear(covariance(a, a), a.variance(a.mean()), 1E-16));
    return 0;
}
