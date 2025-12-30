/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using MatrixType = DenseMatrix<Diff<ScalarType, DiffMode::Reverse, 1>>;
using RandomSource = Random<>;

int main() {
    MatrixType m = MatrixType::random_uniform<RandomSource>(4, 4);
    m.sum().reverse();

    auto v = m.flatten();
    for (size_t i = 0; i < v.getLength(); ++i)
        expect(v.calc(i).grad() == ScalarType(1));
    return 0;
}
