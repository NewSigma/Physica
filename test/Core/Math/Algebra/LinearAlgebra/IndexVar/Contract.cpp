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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;

namespace {
    void tensor_contract() {
        constexpr Var i, j, k, m, n;
        auto B = DenseTensor<T, 2, 2, 2>::random_uniform<RandomSource>(2UZ, 2UZ, 2UZ);
        auto C = DenseTensor<T, 2, 2, 2>::random_uniform<RandomSource>(2UZ, 2UZ, 2UZ);
        DenseTensor<T, 2, 2, 2, 2> A(2UZ, 2UZ, 2UZ, 2UZ);

        A[i, j, m, n] = B[i, j, k] * C[k, m, n];
        for (size_t a = 0; a < 2; ++a)
            for (size_t b = 0; b < 2; ++b)
                for (size_t c = 0; c < 2; ++c)
                    for (size_t d = 0; d < 2; ++d) {
                        T ref = T(0);
                        for (size_t e = 0; e < 2; ++e)
                            ref = fma(B[a, b, e], C[e, c, d], ref);
                        expect(A[a, b, c, d] == ref);
                    }
    }

    void tensor_matrix() {
        constexpr Var i, j, k, m;
        auto tensorB = DenseTensor<T, 2, 2, 2>::random_uniform<RandomSource>(2UZ, 2UZ, 2UZ);
        auto tensorC = DenseMatrix<T>::random_uniform<RandomSource>(2UZ, 2UZ);
        DenseTensor<T, 2, 2, 2> tensorA(2UZ, 2UZ, 2UZ);

        tensorA[i, j, m] = tensorB[i, j, k] * tensorC[k, m];
        for (size_t a = 0; a < 2; ++a)
            for (size_t b = 0; b < 2; ++b)
                for (size_t c = 0; c < 2; ++c) {
                    T ref = T(0);
                    for (size_t d = 0; d < 2; ++d)
                        ref = fma(tensorB[a, b, d], tensorC.calc(d, c), ref);
                    expect(tensorA[a, b, c] == ref);
                }
    }

    void matrix_contract() {
        constexpr Var i, k, m;
        auto C = DenseMatrix<T>::random_uniform<RandomSource>(2UZ, 2UZ);
        DenseMatrix<T> D(2UZ, 2UZ);

        D[i, m] = C[i, k] * C[k, m];
        for (size_t a = 0; a < 2; ++a)
            for (size_t b = 0; b < 2; ++b) {
                T ref = T(0);
                for (size_t d = 0; d < 2; ++d)
                    ref = fma(C.calc(a, d), C.calc(d, b), ref);
                expect(D.calc(a, b) == ref);
            }
    }

    void matrix_vector() {
        constexpr Var i, j;
        auto A = DenseMatrix<T>::random_uniform<RandomSource>(3UZ, 2UZ);
        auto x = DenseVector<T>::random_uniform<RandomSource>(2UZ);
        DenseVector<T> y(3);

        y[i] = A[i, j] * x[j];
        for (size_t r = 0; r < 3; ++r) {
            T ref = T(0);
            for (size_t c = 0; c < 2; ++c)
                ref = fma(A.calc(r, c), x.calc(c), ref);
            expect(y[r] == ref);
        }

        auto ATranspose = DenseMatrix<T>::random_uniform<RandomSource>(2UZ, 3UZ);
        DenseVector<T> y2(3);
        y2[i] = x[j] * ATranspose[j, i];
        for (size_t r = 0; r < 3; ++r) {
            T ref = T(0);
            for (size_t c = 0; c < 2; ++c)
                ref = fma(x.calc(c), ATranspose.calc(c, r), ref);
            expect(y2[r] == ref);
        }
    }

    void vector_outer() {
        constexpr Var i, j;
        auto u = DenseVector<T>::random_uniform<RandomSource>(3UZ);
        auto v = DenseVector<T>::random_uniform<RandomSource>(3UZ);
        DenseMatrix<T> M(3, 3);

        M[i, j] = u[i] * v[j];
        for (size_t a = 0; a < 3; ++a)
            for (size_t b = 0; b < 3; ++b)
                expect(M[a, b] == u[a] * v[b]);
    }

    void tensor_vector() {
        constexpr Var i, j, k, l;
        auto B = DenseTensor<T, 2, 2, 2, 2>::random_uniform<RandomSource>(2UZ, 2UZ, 2UZ, 2UZ);
        auto x = DenseVector<T>::random_uniform<RandomSource>(2UZ);
        DenseTensor<T, 2, 2, 2> A(2UZ, 2UZ, 2UZ);

        A[i, j, k] = B[i, j, k, l] * x[l];
        for (size_t a = 0; a < 2; ++a)
            for (size_t b = 0; b < 2; ++b)
                for (size_t c = 0; c < 2; ++c) {
                    T ref = T(0);
                    for (size_t d = 0; d < 2; ++d)
                        ref = fma(B[a, b, c, d], x[d], ref);
                    expect(A[a, b, c] == ref);
                }
    }

    void tensor_outer() {
        constexpr Var i, j, k, m, n, l;
        auto U = DenseTensor<T, 2, 2, 2>::random_uniform<RandomSource>(2UZ, 2UZ, 2UZ);
        auto V = DenseTensor<T, 2, 2, 2>::random_uniform<RandomSource>(2UZ, 2UZ, 2UZ);
        DenseTensor<T, 2, 2, 2, 2, 2, 2> P(2UZ, 2UZ, 2UZ, 2UZ, 2UZ, 2UZ);

        P[i, j, k, m, n, l] = U[i, j, k] * V[m, n, l];
        for (size_t a = 0; a < 2; ++a)
            for (size_t b = 0; b < 2; ++b)
                for (size_t c = 0; c < 2; ++c)
                    for (size_t d = 0; d < 2; ++d)
                        for (size_t e = 0; e < 2; ++e)
                            for (size_t f = 0; f < 2; ++f)
                                expect(P[a, b, c, d, e, f] == U[a, b, c] * V[d, e, f]);
    }
}

int main() {
    tensor_contract();
    tensor_matrix();
    matrix_contract();
    matrix_vector();
    vector_outer();
    tensor_vector();
    tensor_outer();
    return 0;
}
