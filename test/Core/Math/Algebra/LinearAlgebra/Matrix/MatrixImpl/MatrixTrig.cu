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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.cuh"
#include "Test.h"

using namespace Physica;
using T = float32;
using RandomSource = Random<MCG>;
using M = device_obj<DenseMatrix<T, MatrixMajor::Col>>;

namespace {
    void inverse() noexcept {
        constexpr double Prec = 1E-4;
        const M origin = device_obj<IdentityMatrix<T, 4>>{} + M::random_uniform<RandomSource>(4, 4);
        M inv = origin.triu().inv();
        M prod = origin.triu() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));

        inv = origin.tril().inv();
        prod = origin.tril() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));

        inv = origin.triu_unit().inv();
        prod = origin.triu_unit() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));

        inv = origin.tril_unit().inv();
        prod = origin.tril_unit() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));
    }

    void invGEMM() noexcept {
        constexpr double Prec = 1E-4;
        using Matrix = device_obj<DenseMatrix<T, MatrixMajor::Col>>;
        const Matrix origin = device_obj<IdentityMatrix<T, 4>>{} + Matrix::random_uniform<RandomSource>(4, 4);
        const Matrix rhs = Matrix::random_uniform<RandomSource>(4, 4);
        const Matrix lhs = Matrix::random_uniform<RandomSource>(4, 4);
        Matrix sol(4, 4);
        Matrix prod(4, 4);

        sol = origin.tril().inv() * rhs;
        prod = origin.tril() * sol;
        expect<RandomSource>(matrixNear(prod.toHost(), rhs.toHost(), Prec));

        sol = origin.triu().inv() * rhs;
        prod = origin.triu() * sol;
        expect<RandomSource>(matrixNear(prod.toHost(), rhs.toHost(), Prec));

        sol = origin.tril_unit().inv() * rhs;
        prod = origin.tril_unit() * sol;
        expect<RandomSource>(matrixNear(prod.toHost(), rhs.toHost(), Prec));

        sol = origin.triu_unit().inv() * rhs;
        prod = origin.triu_unit() * sol;
        expect<RandomSource>(matrixNear(prod.toHost(), rhs.toHost(), Prec));

        sol = lhs * origin.tril().inv();
        prod = sol * origin.tril();
        expect<RandomSource>(matrixNear(prod.toHost(), lhs.toHost(), Prec));

        sol = lhs * origin.triu().inv();
        prod = sol * origin.triu();
        expect<RandomSource>(matrixNear(prod.toHost(), lhs.toHost(), Prec));

        sol = lhs * origin.tril_unit().inv();
        prod = sol * origin.tril_unit();
        expect<RandomSource>(matrixNear(prod.toHost(), lhs.toHost(), Prec));

        sol = lhs * origin.triu_unit().inv();
        prod = sol * origin.triu_unit();
        expect<RandomSource>(matrixNear(prod.toHost(), lhs.toHost(), Prec));
    }
}

int main() {
    inverse();
    invGEMM();
    return 0;
}
