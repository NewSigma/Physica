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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MCG>;
using T = float64;

int main() {
    const auto a = MatrixND<T>::random_uniform<RandomSource>(4, 4);
    const auto b = MatrixND<T>::random_uniform<RandomSource>(6, 6);
    const auto v = VectorND<T>::random_uniform<RandomSource>(kronecker(a, b).getRow());
    expect<RandomSource>(vectorNear(MatrixND<T>(kronecker(a, b)) * v, kronecker(a, b) * v, 3UL));
    return 0;
}
