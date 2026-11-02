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
    void testSize() noexcept {
        // Test that getRow() and getCol() are correct for rectangle matrix
        M rectangle(5, 4);
        expect(rectangle.triu().getRow() == 4);
        expect(rectangle.triu().getCol() == 4);
        expect(rectangle.transpose().triu().getRow() == 4);
        expect(rectangle.transpose().triu().getCol() == 5);

        expect(rectangle.tril().getRow() == 5);
        expect(rectangle.tril().getCol() == 4);
        expect(rectangle.transpose().tril().getRow() == 4);
        expect(rectangle.transpose().tril().getCol() == 4);
    }

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
        const auto check = [&](const Matrix auto& trig) {
            M inv = trig.inv();
            M prod = trig * inv;
            expect<RandomSource>(matrixNear(prod.toHost(), IdentityMatrix<T, 4>(4), Prec));
        };
        check(data.triu());
        check(data.tril());
        check(data.triu_unit());
        check(data.tril_unit());
        check(data.transpose().triu());
        check(data.transpose().tril());
        check(data.transpose().triu_unit());
        check(data.transpose().tril_unit());
    }

    void invGEMV() noexcept {
        constexpr double Prec = 1E-4;
        using DVector = device_obj<VectorND<T>>;
        const M data = device_obj<IdentityMatrix<T, 4>>{} + M::random_uniform<RandomSource>(4, 4);
        const VectorND<T> hostV = VectorND<T>::random_uniform<RandomSource>(4) + T(0.5); // Add 0.5 so a reference near zero inflates the error
        const DVector v = hostV.toDevice();
        DVector sol;
        const auto check = [&](const Matrix auto& trig) {
            sol = trig.inv() * v;
            const DVector prod = trig * sol;
            expect<RandomSource>(vectorNear(prod.toHost(), v.toHost(), Prec));
        };

        check(data.tril());
        check(data.tril_unit());
        check(data.triu());
        check(data.triu_unit());
        check(data.transpose().tril());
        check(data.transpose().triu());
        check(data.transpose().tril_unit());
        check(data.transpose().triu_unit());
    }

    void invGEMM() noexcept {
        constexpr double Prec = 1E-4;
        const MatrixND<T> hostRhs = MatrixND<T>::random_uniform<RandomSource>(4, 4) + T(0.5);
        const MatrixND<T> hostLhs = MatrixND<T>::random_uniform<RandomSource>(4, 4) + T(0.5);
        const M rhs = hostRhs.toDevice();
        const M lhs = hostLhs.toDevice();
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
    testSize();
    gemm();
    inverse();
    invGEMV();
    invGEMM();
    return 0;
}
