/*
 * Copyright 2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
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

    auto samples = VectorND<T>::generate(NumSystem, [&](int) {
        return idet.compute_base<RandomSource>(A, 0.05, NumSample);
    });
    auto mean = samples.mean();
    auto devia = samples.deviation(mean) / sqrt(T(samples.getLength()));
    expect<RandomSource>(abs(mean - A.lnAbsDet()) < T(3) * devia);
    expect<RandomSource>(devia < abs(mean) * T(1E-2));
    return 0;
}
