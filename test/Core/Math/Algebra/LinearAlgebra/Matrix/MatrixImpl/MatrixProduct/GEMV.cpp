/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Scalar/Diff.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void forward() {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Forward>;
        auto m = MatrixND<dfloat>::random_uniform<RandomSource>(4, 4);
        auto x = VectorND<dfloat>::random_uniform<RandomSource>(4);
        auto y = VectorND<dfloat>::junk(4);
        y = m * x;
        expect<RandomSource>(vectorNear(y.values(), VectorND<T>(m.values() * x.values()), 1E-13));
        expect(y.grads().isZero());
    }

    void reverse() {
        using ScalarType = Diff<float32, DiffMode::Reverse>;
        using MatrixType = DenseMatrix<float32, MatrixMajor::Col, 3, 3>;
        using DVector = Vector3D<ScalarType>;
        using DMatrix = DenseMatrix<ScalarType, MatrixMajor::Col, 3, 3>;
        DMatrix m{1, 2, 3, 4, 5, 6, 7, 8, 9};
        DVector x{1, 2, 3};
        {
            CoDiff<DVector> y = m * x;
            y.sum().reverse();
        }
        const Vector3D<float32> dx = x.grads();
        for (size_t i = 0; i < dx.getLength(); ++i)
            expect(scalarNear(dx[i], m.values().col(i).sum(), 1E-15));

        const MatrixType dm = m.grads();
        for (size_t r = 0; r < dm.getRow(); ++r)
            for (size_t c = 0; c < dm.getCol(); ++c)
                expect(scalarNear(dm[r, c], x.values().calc(c), 1E-15));
    }

    void transpose() {
        using T = float64;
        constexpr size_t Row = 64;
        constexpr size_t Col = 32;
        auto check = []<Matrix M>(const M& mat) {
            const size_t row = mat.getRow();
            const size_t col = mat.getCol();
            const auto v = VectorND<T>::random_uniform<RandomSource>(col);
            VectorND<T> expected(row);
            for (size_t i = 0; i < row; ++i) {
                T sum = 0;
                for (size_t j = 0; j < col; ++j)
                    sum += mat.calc(i, j) * v[j];
                expected[i] = sum;
            }
            VectorND<T> actual = mat * v;
            expect<RandomSource>(vectorNear(actual, expected, 1E-12));
        };

        auto rowMajor = DenseMatrix<T, MatrixMajor::Row>::random_uniform<RandomSource>(Row, 2 * Col);
        auto colMajor = DenseMatrix<T, MatrixMajor::Col>::random_uniform<RandomSource>(Row, 2 * Col);
        /* Whole matrix */
        check(rowMajor);
        check(rowMajor.transpose());
        check(colMajor);
        check(colMajor.transpose());
        /* Strided block, as is frequently produced by column slicing */
        check(rowMajor.leftCols(Col));
        check(rowMajor.leftCols(Col).transpose());
        check(colMajor.leftCols(Col));
        check(colMajor.leftCols(Col).transpose());
    }
}

int main() {
    forward();
    reverse();
    transpose();
    return 0;
}
