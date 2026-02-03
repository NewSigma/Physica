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

int main() {
    using T = float32;
    using Matrix4D = DenseMatrix<T, MatrixMajor::Col, 4, 4>;
    auto continuous = Matrix4D::random_uniform<Random<>>(4, 4);
    auto lMatrix = continuous.topLeftCorner(3);
    for (int i = 0; i < 3; ++i)
        expect(lMatrix.diag(0)[i] == continuous.diag()[i]); // 0 sub-diagonal is essentially diagonal
    return 0;
}
