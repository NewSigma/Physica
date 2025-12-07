/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/LUDecomp.h"

using namespace Physica;
using Matrix3D = DenseMatrix<float64, MatrixOption::Col, 3, 3>;
using Matrix4D = DenseMatrix<float64, MatrixOption::Col, 4, 4>;

namespace {
    template<bool Pivot>
    void test(const Matrix auto& answer, double prec) {
        using M = std::remove_cvref<decltype(answer)>::type;
        LUDecomp<float64, Pivot> lu(answer.getRow());
        auto product = [&, prec]() {
            M matrixL = lu.getMatrixLU().tril();
            matrixL.diag() = float64(1);
            M result = matrixL * lu.getMatrixU();
            if constexpr (Pivot)
                result = M(lu.getPerm() * result);

            if (!matrixNear(result, answer, prec))
                exit(1);
        };

        lu.compute_base(answer);
        product();

        if constexpr (HasMKL()) {
            lu.compute_mkl(answer);
            product();
        }
    }
}

int main() {
    const Matrix3D mat1{{2, 3, 4}, {1, 1, 9}, {1, 2, -6}};
    test<false>(mat1, 1E-15);
    test<true>(mat1, 1E-15);
    test<true>(Matrix4D{
        {1,  2, 0, 1},
        {0,  1, 1, 0},
        {2,  0, 1, 1},
        {1,  1, 0, 1}
    }, 1E-15);
    return 0;
}
