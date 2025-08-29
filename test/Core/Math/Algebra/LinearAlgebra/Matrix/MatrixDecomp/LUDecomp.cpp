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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/LUDecomp.h"

using namespace Physica;
using Matrix3D = DenseMatrix<float64, MatrixOption::Col | MatrixOption::Vector, 3, 3>;

void test(const Matrix auto& answer) {
    LUDecomp<float64, false> lu(answer.getRow());
    auto product = [&]() {
        Matrix3D matrixL = lu.getMatrixLU().tril();
        matrixL.diag() = float64(1);
        Matrix3D result = matrixL * lu.getMatrixU();
        if (!matrixNear(result, answer, 1E-15))
            exit(1);
    };

    lu.compute_base(answer);
    product();

    if constexpr (HasMKL()) {
        lu.compute_mkl(answer);
        product();
    }
}

int main() {
    const Matrix3D mat1{{2, 3, 4}, {1, 1, 9}, {1, 2, -6}};
    test(mat1);
    return 0;
}
