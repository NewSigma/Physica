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
using RandomSource = Random<>;

namespace {
    void sum() {
        constexpr int N = 16;
        using T = float64;
        const auto a = VectorND<T>::random_uniform<RandomSource>(N);
        const auto b = VectorND<T>::random_uniform<RandomSource>(N);
        const MatrixND<T> y = a * b.transpose();
        expect<RandomSource>(scalarNear((a * b.transpose()).sum(), y.sum(), 1E-15));
    }
}

int main() {
    sum();
    return 0;
}
