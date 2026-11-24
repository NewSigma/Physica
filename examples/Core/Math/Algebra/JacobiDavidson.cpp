/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using MatrixType = DenseHermiteMatrix<cfloat64>;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

int main() {
    const auto mat = MatrixType::random_uniform<RandomSource>(5000);

    JacobiDavidson<cfloat64> jd(mat.getRow(), 4);
    jd.compute(mat, VectorND<cfloat64>::random_uniform<RandomSource>(mat.getRow()));
    jd.sort();

    std::cout << jd.getEigenvalues().reals().format() << '\n';
    return 0;
}
