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
using RandomSource = Random<MCG>;
using Matrix4D = DenseMatrix<T, MatrixMajor::Col, 4, 4>;

namespace {
    void testSize() noexcept {
        // Test that getRow() and getCol() are correct for rectangle matrix
        MatrixND<T> rectangle(5, 4);
        expect(rectangle.triu().getRow() == 4);
        expect(rectangle.triu().getCol() == 4);
        expect(rectangle.transpose().triu().getRow() == 4);
        expect(rectangle.transpose().triu().getCol() == 5);

        expect(rectangle.tril().getRow() == 5);
        expect(rectangle.tril().getCol() == 4);
        expect(rectangle.transpose().tril().getRow() == 4);
        expect(rectangle.transpose().tril().getCol() == 4);
    }

    void invAndProd() noexcept {
        constexpr double Prec = 1E-12;
        auto origin = Matrix4D::random_uniform<RandomSource>(4, 4);
        Matrix4D inv = origin.triu().inv();
        Matrix4D prod = origin.triu() * inv;
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));
        prod = inv * origin.triu();
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));

        inv = origin.tril().inv();
        prod = origin.tril() * inv;
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));
        prod = inv * origin.tril();
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));

        inv = origin.triu_unit().inv();
        prod = origin.triu_unit() * inv;
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));
        prod = inv * origin.triu_unit();
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));

        inv = origin.tril_unit().inv();
        prod = origin.tril_unit() * inv;
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));
        prod = inv * origin.tril_unit();
        expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));
    }

    void invGEMV() {
        constexpr double Prec = 1E-12;
        auto m = Matrix4D::random_uniform<RandomSource>(4, 4);
        auto v = Vector4D<T>::random_uniform<RandomSource>(4);
        Vector4D<T> sol = m.tril().inv() * v;
        expect<RandomSource>(vectorNear(m.tril() * sol, v, Prec));

        sol = m.tril_unit().inv() * v;
        expect<RandomSource>(vectorNear(m.tril_unit() * sol, v, Prec));

        sol = m.triu().inv() * v;
        expect<RandomSource>(vectorNear(m.triu() * sol, v, Prec));

        sol = m.triu_unit().inv() * v;
        expect<RandomSource>(vectorNear(m.triu_unit() * sol, v, Prec));
    }
}

int main() {
    testSize();
    invAndProd();
    invGEMV();
    return 0;
}
