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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.h"
#include "Test.h"

using namespace Physica;

namespace {
    void trivial_inv() noexcept {
        using T = float64;
        using Matrix4D = DenseMatrix<T, MatrixOption::Col, 4, 4>;
        auto origin = Matrix4D::random_uniform<Random<MCG, 1234>>(4, 4);
        Matrix4D inv = origin.triu().inv();
        Matrix4D prod = origin.triu() * inv;
        expect(matrixNear(prod, IdentityMatrix<T, 4>(4), 1E-14));
    }
}

int main() {
    trivial_inv();
    return 0;
}
