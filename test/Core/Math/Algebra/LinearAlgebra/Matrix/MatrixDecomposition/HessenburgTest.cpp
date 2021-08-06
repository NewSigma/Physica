/*
 * Copyright 2021 WeiBo He.
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
#include "TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Hessenburg.h"
#include "iostream"
using namespace Physica::Core;

int main() {
    using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 4, 4>;
    const MatrixType mat{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    const MatrixType answer{{1, -5.38516, 0, 0},
                            {-16.5269, 33.8276, -6.75721, 0},
                            {-1.36458, 0.591256, -0.827586, 0},
                            {0, 0, 0, 0}};
    MatrixType result = Hessenburg(mat);
    std::cout << result << std::endl;
    std::cout << answer << std::endl;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (!floatNear(result(i, j), answer(i, j), 1E-5)) {
                std::cout << result(i, j) << ' ' << answer(i, j) << '\n';
                return 1;
            }
    return 0;
}
