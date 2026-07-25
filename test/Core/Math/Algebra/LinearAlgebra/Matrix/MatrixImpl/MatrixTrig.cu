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
    void gemm() noexcept {
        constexpr double Prec = 1E-4;
        const M rhs = M::random_uniform<RandomSource>(4, 4);
        const M lhs = M::random_uniform<RandomSource>(4, 4);
        M result, answer, dense;
        const auto doubleProd = [&](const Matrix auto& trig) {
            dense = trig;
            result = trig * rhs;
            answer = dense * rhs;
            expect<RandomSource>(matrixNear(result.toHost(), answer.toHost(), Prec));

            result = lhs * trig;
            answer = lhs * dense;
            expect<RandomSource>(matrixNear(result.toHost(), answer.toHost(), Prec));
        };

        const M data = M::random_uniform<RandomSource>(4, 4);
        doubleProd(data.tril());
        doubleProd(data.triu());
        doubleProd(data.tril_unit());
        doubleProd(data.triu_unit());
        doubleProd(data.transpose().tril());
        doubleProd(data.transpose().triu());
        doubleProd(data.transpose().tril_unit());
        doubleProd(data.transpose().triu_unit());
    }

    void inverse() noexcept {
        constexpr double Prec = 1E-4;
        const M data = device_obj<IdentityMatrix<T, 4>>{} + M::random_uniform<RandomSource>(4, 4);
        M inv = data.triu().inv();
        M prod = data.triu() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));

        inv = data.tril().inv();
        prod = data.tril() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));

        inv = data.triu_unit().inv();
        prod = data.triu_unit() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));

        inv = data.tril_unit().inv();
        prod = data.tril_unit() * inv;
        expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));
    }

    void invGEMM() noexcept {
        constexpr double Prec = 1E-4;
        const M rhs = M::random_uniform<RandomSource>(4, 4);
        const M lhs = M::random_uniform<RandomSource>(4, 4);
        M sol(4, 4);
        M prod(4, 4);
        const auto doubleProd = [&](const Matrix auto& trig) {
            sol = trig.inv() * rhs;
            prod = trig * sol;
            expect<RandomSource>(matrixNear(prod.toHost(), rhs.toHost(), Prec));

            sol = lhs * trig.inv();
            prod = sol * trig;
            expect<RandomSource>(matrixNear(prod.toHost(), lhs.toHost(), Prec));
        };

        const M data = device_obj<IdentityMatrix<T, 4>>{} + M::random_uniform<RandomSource>(4, 4);
        doubleProd(data.tril());
        doubleProd(data.triu());
        doubleProd(data.tril_unit());
        doubleProd(data.triu_unit());
        doubleProd(data.transpose().tril());
        doubleProd(data.transpose().triu());
        doubleProd(data.transpose().tril_unit());
        doubleProd(data.transpose().triu_unit());
    }
}

int main() {
    gemm();
    inverse();
    invGEMM();
    return 0;
}
