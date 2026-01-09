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
using Matrix3D = DenseMatrix<float64, MatrixOption::Col, 3, 3>;
using Matrix4D = DenseMatrix<float64, MatrixOption::Col, 4, 4>;

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
}

int main() {
    const Matrix3D mat1{{2, 3, 4}, {1, 1, 9}, {1, 2, -6}};
    const Matrix4D mat2{
        {1,  2, 0, 1},
        {0,  1, 1, 0},
        {2,  0, 1, 1},
        {1,  1, 0, 1}
    };
    decomp<false>(mat1, 1E-15);
    decomp<true>(mat1, 1E-15);
    decomp<true>(mat2, 1E-15);

    inverse<false>(mat1, 1E-15);
    inverse<true>(mat1, 1E-14);
    inverse<true>(mat2, 1E-10);
    return 0;
}
