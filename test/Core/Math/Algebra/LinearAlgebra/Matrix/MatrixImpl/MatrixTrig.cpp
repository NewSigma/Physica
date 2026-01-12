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
using T = float64;
using RandomSource = Random<MCG, 1234>;
using Matrix4D = DenseMatrix<T, MatrixOption::Col, 4, 4>;

namespace {
    void trivial_inv() noexcept {
        auto origin = Matrix4D::random_uniform<RandomSource>(4, 4);
        Matrix4D inv = origin.triu().inv();
        Matrix4D prod = origin.triu() * inv;
        expect(matrixNear(prod, IdentityMatrix<T, 4>(4), 1E-14));

        inv = origin.tril().inv();
        prod = origin.tril() * inv;
        expect(matrixNear(prod, IdentityMatrix<T, 4>(4), 1E-14));

        inv = origin.triu_unit().inv();
        prod = origin.triu_unit() * inv;
        expect(matrixNear(prod, IdentityMatrix<T, 4>(4), 1E-14));

        inv = origin.tril_unit().inv();
        prod = origin.tril_unit() * inv;
        expect(matrixNear(prod, IdentityMatrix<T, 4>(4), 1E-14));
    }

    void invGEMV() {
        auto m = Matrix4D::random_uniform<RandomSource>(4, 4);
        auto v = Vector4D<T>::random_uniform<RandomSource>(4);
        Vector4D<T> sol = m.tril().inv() * v;
        expect(vectorNear(m.tril() * sol, v, 1E-15));

        sol = m.tril_unit().inv() * v;
        expect(vectorNear(m.tril_unit() * sol, v, 1E-15));
    }
}

int main() {
    trivial_inv();
    invGEMV();
    return 0;
}
