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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;

int main() {
    auto matD = DiagMatrix<T>(VectorND<T>::random_uniform<Random<>>(16));
    auto x = VectorND<T>::random_uniform<Random<>>(16);
    VectorND<T> result = hadamard(x, matD.diag());
    x = matD * x;
    expect(x == result); // Test that there is no aliasing problem
    return 0;
}
