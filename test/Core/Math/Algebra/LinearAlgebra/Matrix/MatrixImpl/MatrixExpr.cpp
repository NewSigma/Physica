/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;

namespace {
    void general() {
        DenseMatrix<float64, MatrixMajor::Col, 3, 3> mat1{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
        DenseMatrix<float32, MatrixMajor::Col, 3, 3> mat2{1, 1, 1, 1, 1, 1, 1, 1, 1};
        {
            auto mat = DenseMatrix<float64, MatrixMajor::Row, 3, 3>::junk(3, 3);
            mat = -(mat1 + mat2);
            for (size_t i = 0; i < mat.getRow(); ++i)
                for (size_t j = 0; j < mat.getCol(); ++j)
                    expect(mat[i, j] == float64(-2));
        }
        {
            auto mat = DenseMatrix<float64, MatrixMajor::Row, 3, 3>::junk(3, 3);
            mat = mat1 * mat2;
            for (size_t i = 0; i < mat.getRow(); ++i)
                for (size_t j = 0; j < mat.getCol(); ++j)
                    expect(mat[i, j] == float64(3));
        }
    }

    void block1x1() {
        using ScalarType = float64;
        DenseMatrix<ScalarType, MatrixMajor::Col, Physica::Dynamic, 1> mat(2, 1);
        mat[0, 0] = 1.0;
        mat[1, 0] = 2.0;
        expect(mat.row(1)[0] == ScalarType(2));
    }

    void lifetimeCXX23() {
        // Test that we are free of a set of lifetime problems under CXX23
        using T = float32;
        auto diag = (MatrixND<T>::identity(3) + MatrixND<T>(3, 3, 5)).diag();
        for (int i = 0; i < 3; ++i)
            expect(diag.calc(i) == T(6));
    }

    void simplify() {
        syntax_only([]() {
            using T = float32;
            auto x = -(MatrixND<T>() * T(2));
            static_assert(x.getExprID() == ExprID::Mul);
        });
    }

    void hermite() {
        // Elementwise regression test
        using T = cfloat64;
        using MatrixType = DenseMatrix<T, MatrixMajor::Col, 2, 2>;
        const MatrixType mat1{{T(1, 1), T(2, 3)}, {T(4, 5), T(6, 7)}};
        const MatrixType mat2{{T(8, 9), T(10, 11)}, {T(12, 13), T(14, 15)}};
        const MatrixType result = (mat1 + mat2).hermite();
        for (size_t i = 0; i < result.getRow(); ++i)
            for (size_t j = 0; j < result.getCol(); ++j)
                expect(result[i, j] == (mat1.calc(j, i) + mat2.calc(j, i)).conjugate());
    }
}

int main() {
    general();
    block1x1();
    lifetimeCXX23();
    simplify();
    hermite();
    return 0;
}
