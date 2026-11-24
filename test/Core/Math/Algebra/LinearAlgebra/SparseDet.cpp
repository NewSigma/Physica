/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/SparseDet.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;
constexpr int NumSample = 512;
constexpr int NumSystem = 256;

int main() {
    const DenseMatrix<T, MatrixMajor::Row, 3, 3> A{
        {3,  2, 0},
        {1, -1, 0},
        {0,  5, 1}
    };
    SparseDet<T> idet(GMRES<T>(3, 3));
    idet.getGMRES().setTolerance(1E-4);

    auto samples = VectorND<T>::generate([&](int) {
        return idet.compute_base<RandomSource>(A, 0.05, NumSample);
    }, NumSystem);
    auto mean = samples.mean();
    auto devia = samples.deviation(mean) / sqrt(T(samples.getLength()));
    expect<RandomSource>(abs(mean - A.lnAbsDet()) < T(3) * devia);
    expect<RandomSource>(devia < abs(mean) * T(1E-2));
    return 0;
}
