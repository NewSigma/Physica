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
