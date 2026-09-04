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

    void gemm() {
        constexpr double Prec = 1E-11;
        const auto rhs = Matrix4D::random_uniform<RandomSource>(4, 4);
        const auto lhs = Matrix4D::random_uniform<RandomSource>(4, 4);
        Matrix4D dense, result, answer;
        const auto doubleProd = [&](const Matrix auto& trig) {
            dense = trig;
            answer = dense * rhs;
            (trig * rhs).assign_base(result);
            expect<RandomSource>(matrixNear(result, answer, Prec));
            if constexpr (HasMKL()) {
                (trig * rhs).assign_mkl(result);
                expect<RandomSource>(matrixNear(result, answer, Prec));
            }

            answer = lhs * dense;
            (lhs * trig).assign_base(result);
            expect<RandomSource>(matrixNear(result, answer, Prec));
            if constexpr (HasMKL()) {
                (lhs * trig).assign_mkl(result);
                expect<RandomSource>(matrixNear(result, answer, Prec));
            }
        };

        const auto data = Matrix4D::random_uniform<RandomSource>(4, 4);
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
        constexpr double Prec = 1E-10;
        const auto data = Matrix4D::template random_uniform<RandomSource>(4, 4);
        const auto check = [&](const Matrix auto& trig) {
            Matrix4D inv = trig.inv();
            Matrix4D prod = trig * inv;
            expect<RandomSource>(matrixNear(prod, IdentityMatrix<T, 4>(4), Prec));
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

    void invGEMM() {
        constexpr double Prec = 1E-11;
        const auto rhs = Matrix4D::random_uniform<RandomSource>(4, 4);
        const auto lhs = Matrix4D::random_uniform<RandomSource>(4, 4);
        Matrix4D sol, prod;
        const auto doubleProd = [&](const Matrix auto& trig) {
            Matrix4D sol = trig.inv() * rhs;
            Matrix4D prod = trig * sol;
            expect<RandomSource>(matrixNear(prod, rhs, Prec));

            sol = lhs * trig.inv();
            prod = sol * trig;
            expect<RandomSource>(matrixNear(prod, lhs, Prec));
        };

        const auto data = Matrix4D::random_uniform<RandomSource>(4, 4);
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
