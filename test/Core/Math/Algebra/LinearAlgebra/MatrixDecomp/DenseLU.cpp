/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.h"
#include "Test.h"

using namespace Physica;
using Matrix3D = DenseMatrix<float64, MatrixMajor::Col, 3, 3>;
using Matrix4D = DenseMatrix<float64, MatrixMajor::Col, 4, 4>;

namespace {
    template<bool Pivot>
    void decomp(const Matrix auto& answer, double prec) {
        using M = std::remove_cvref<decltype(answer)>::type;
        DenseLU<float64, Pivot> lu(answer.getRow());
        auto product = [&, prec]() {
            M matrixL = lu.getMatrixLU().tril();
            matrixL.diag() = float64(1);
            M result = matrixL * lu.getMatrixU();
            if constexpr (Pivot)
                result = M(lu.getPerm() * result);
            expect(matrixNear(result, answer, prec));
        };

        lu.compute_base(answer);
        product();

        if constexpr (HasMKL()) {
            lu.compute_mkl(answer);
            product();
        }
    }

    template<bool Pivot>
    void inverse(const Matrix auto& source, double prec) {
        using M = std::remove_cvref<decltype(source)>::type;
        using T = M::ScalarType;
        DenseLU<T, Pivot> lu(source.getRow());
        lu.compute(source);

        M buffer, prod;
        buffer.resize(source);
        if constexpr (!Pivot) {
            lu.inv().assign_base(buffer);
            prod = buffer * source;
            expect(matrixNear(prod, IdentityMatrix<T>(source.getRow()), prec));
        }

        if constexpr (HasMKL()) {
            lu.inv().assign_mkl(buffer);
            prod = buffer * source;
            expect(matrixNear(prod, IdentityMatrix<T>(source.getRow()), prec));
        }
    }

    template<bool Pivot>
    void inverseGEMV(double prec) {
        using T = float64;
        using Matrix4D = DenseMatrix<T, MatrixMajor::Row, 4, 4>;
        const Matrix4D A{
                {-0.000696013585639699,     0.816492585748236, 0.0216969440126965, -0.0884307621566726},
                {    0.691809621910274, -0.000696013585639699,  0.131671000379563, -0.0701048797366553},
                {  -0.0701048797366553,   -0.0884307621566726, -0.131016640264434,   0.788769710999288},
                {    0.131671000379563,    0.0216969440126965,  0.643819646681362,  -0.131016640264434}
        };
        const Vector4D<T> b{4.316511702487202E-1, 1.548712563601895E-2, 9.840637243791538E-1, 1.671684099146560E-1};
        const Vector4D<T> answer{0.06910464034803039, 0.6682416388355244, 0.5106380624075890, 1.413471488683768};
        auto lu = DenseLU<T, Pivot>(A);
        const Vector4D<T> result = lu.inv() * b;
        expect(vectorNear(result, answer, prec));
    }
}

int main() {
    {
        // We know this matrix does not need pivoting
        const Matrix3D mat{
                {2, 3,  4},
                {1, 1,  9},
                {1, 2, -6}
        };
        decomp<false>(mat, 1E-15);
        inverse<false>(mat, 1E-15);
    }
    {
        // Otherwise, always prefer pivoting
        const Matrix4D mat1{
                {1, 2, 0, 1},
                {0, 1, 1, 0},
                {2, 0, 1, 1},
                {1, 1, 0, 1}
        };
        decomp<true>(mat1, 1E-15);
        inverse<true>(mat1, 1E-14);

        // Partial pivoting does not work for it
        const Matrix4D mat2{
                {1,  1,  1,  1},
                {1,  1, -1, -1},
                {1, -1,  1, -1},
                {1, -1, -1,  1}
        };
        decomp<true>(mat2, 1E-15);
        inverse<true>(mat2, 1E-10);
    }
    inverseGEMV<false>(1E-11); // Precision is lower without pivoting
    inverseGEMV<true>(1E-15);
    return 0;
}
