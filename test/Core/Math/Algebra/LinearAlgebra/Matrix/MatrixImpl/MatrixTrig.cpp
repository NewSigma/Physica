/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        // Empirically I + U[0,1) keeps cond <= 14, so the fixed tolerance never has to be relaxed.
        constexpr double Prec = 1E-10;
        const Matrix4D data = IdentityMatrix<T, 4>(4) + Matrix4D::template random_uniform<RandomSource>(4, 4);
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
        constexpr double Prec = 1E-10;
        const Matrix4D m = IdentityMatrix<T, 4>(4) + Matrix4D::random_uniform<RandomSource>(4, 4);
        Vector4D<T> v = Vector4D<T>::random_uniform<RandomSource>(4) + T(0.5); // Add 0.5 so a reference near zero inflates the error
        Vector4D<T> sol;
        const auto check = [&](const Matrix auto& trig) {
            (trig.inv() * v).assign_base(sol);
            expect<RandomSource>(vectorNear(trig * sol, v, Prec));

            if constexpr (HasMKL()) {
                (trig.inv() * v).assign_mkl(sol);
                expect<RandomSource>(vectorNear(trig * sol, v, Prec));
            }
        };

        check(m.tril());
        check(m.tril_unit());
        check(m.triu());
        check(m.triu_unit());
        check(m.transpose().tril());
        check(m.transpose().triu());
        check(m.transpose().tril_unit());
        check(m.transpose().triu_unit());
    }

    void invGEMM() {
        constexpr double Prec = 1E-10;
        Matrix4D rhs = Matrix4D::random_uniform<RandomSource>(4, 4) + T(0.5);
        Matrix4D lhs = Matrix4D::random_uniform<RandomSource>(4, 4) + T(0.5);
        Matrix4D sol, prod;
        const auto doubleProd = [&](const Matrix auto& trig) {
            sol = trig.inv() * rhs;
            prod = trig * sol;
            expect<RandomSource>(matrixNear(prod, rhs, Prec));

            sol = lhs * trig.inv();
            prod = sol * trig;
            expect<RandomSource>(matrixNear(prod, lhs, Prec));
        };

        const Matrix4D data = IdentityMatrix<T, 4>(4) + Matrix4D::random_uniform<RandomSource>(4, 4);
        doubleProd(data.tril());
        doubleProd(data.triu());
        doubleProd(data.tril_unit());
        doubleProd(data.triu_unit());
        doubleProd(data.transpose().tril());
        doubleProd(data.transpose().triu());
        doubleProd(data.transpose().tril_unit());
        doubleProd(data.transpose().triu_unit());
    }
    /**
     * A non-square operand tells the row and column counts apart, which a wrong transpose flag or a wrong
     * leading dimension swaps silently. Row-major layouts exercise the other cblas layout.
     */
    template<int DataMajor, int TargetMajor>
    void gemmNonsquare() {
        constexpr size_t Order = 8;
        constexpr size_t Other = 16;
        constexpr double Prec = 1E-11;
        using Data = DenseMatrix<T, DataMajor>;
        using Target = DenseMatrix<T, TargetMajor>;

        const auto data = Data::template random_uniform<RandomSource>(Order, Order);
        const auto rhs = Data::template random_uniform<RandomSource>(Order, Other);
        const auto lhs = Data::template random_uniform<RandomSource>(Other, Order);
        Target dense(Order, Order);
        Target answer(Order, Other);
        Target result(Order, Other);
        Target answer2(Other, Order);
        Target result2(Other, Order);

        const auto check = [&](const Matrix auto& trig) {
            dense = trig;
            answer = dense * rhs;
            (trig * rhs).assign(result);
            expect<RandomSource>(matrixNear(result, answer, Prec));
            (trig * rhs).assign_base(result);
            expect<RandomSource>(matrixNear(result, answer, Prec));

            answer2 = lhs * dense;
            (lhs * trig).assign(result2);
            expect<RandomSource>(matrixNear(result2, answer2, Prec));
            (lhs * trig).assign_base(result2);
            expect<RandomSource>(matrixNear(result2, answer2, Prec));
        };

        check(data.tril());
        check(data.triu());
        check(data.tril_unit());
        check(data.triu_unit());
        check(data.transpose().tril());
        check(data.transpose().triu());
        check(data.transpose().tril_unit());
        check(data.transpose().triu_unit());
    }

    template<int DataMajor, int TargetMajor>
    void invGEMMNonsquare() {
        constexpr size_t Order = 8;
        constexpr size_t Other = 16;
        constexpr double Prec = 1E-10;
        using Data = DenseMatrix<T, DataMajor>;
        using Target = DenseMatrix<T, TargetMajor>;

        auto data = Data::template random_uniform<RandomSource>(Order, Order);
        for (size_t i = 0; i < Order; ++i)
            data[i, i] = data[i, i] + T(4); // Diagonally dominant, so that the inverse stays well conditioned
        const auto rhs = Data::template random_uniform<RandomSource>(Order, Other) + T(0.5);
        const auto lhs = Data::template random_uniform<RandomSource>(Other, Order) + T(0.5);
        Target sol(Order, Other);
        Target prod(Order, Other);
        Target sol2(Other, Order);
        Target prod2(Other, Order);

        const auto check = [&](const Matrix auto& trig) {
            sol = trig.inv() * rhs;
            prod = trig * sol;
            expect<RandomSource>(matrixNear(prod, rhs, Prec));

            sol2 = lhs * trig.inv();
            prod2 = sol2 * trig;
            expect<RandomSource>(matrixNear(prod2, lhs, Prec));
        };

        check(data.tril());
        check(data.triu());
        check(data.tril_unit());
        check(data.triu_unit());
        check(data.transpose().tril());
        check(data.transpose().triu());
        check(data.transpose().tril_unit());
        check(data.transpose().triu_unit());
    }

    template<int DataMajor>
    void invGEMVRowMajor() {
        constexpr size_t Order = 8;
        constexpr double Prec = 1E-10;
        using Data = DenseMatrix<T, DataMajor>;

        auto m = Data::template random_uniform<RandomSource>(Order, Order);
        for (size_t i = 0; i < Order; ++i)
            m[i, i] = m[i, i] + T(4);
        const VectorND<T> v = VectorND<T>::template random_uniform<RandomSource>(Order) + T(0.5);
        VectorND<T> sol(Order);
        const auto check = [&](const Matrix auto& trig) {
            sol = trig.inv() * v;
            expect<RandomSource>(vectorNear(trig * sol, v, Prec));
        };

        check(m.tril());
        check(m.tril_unit());
        check(m.triu());
        check(m.triu_unit());
        check(m.transpose().tril());
        check(m.transpose().triu());
        check(m.transpose().tril_unit());
        check(m.transpose().triu_unit());
    }
}

int main() {
    testSize();
    gemm();
    inverse();
    invGEMV();
    invGEMM();
    gemmNonsquare<MatrixMajor::Col, MatrixMajor::Col>();
    gemmNonsquare<MatrixMajor::Col, MatrixMajor::Row>();
    gemmNonsquare<MatrixMajor::Row, MatrixMajor::Col>();
    gemmNonsquare<MatrixMajor::Row, MatrixMajor::Row>();
    invGEMMNonsquare<MatrixMajor::Col, MatrixMajor::Col>();
    invGEMMNonsquare<MatrixMajor::Col, MatrixMajor::Row>();
    invGEMMNonsquare<MatrixMajor::Row, MatrixMajor::Col>();
    invGEMMNonsquare<MatrixMajor::Row, MatrixMajor::Row>();
    invGEMVRowMajor<MatrixMajor::Row>();
    return 0;
}
