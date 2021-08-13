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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/RealSchur.h"

using namespace Physica::Core;

int main() {
    using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3>;
    const MatrixType mat{{-149, 537, -27}, {-50, 180, -9}, {-154, 546, -25}};
    const MatrixType answer{{1, 0, 0}, {-7.1119, 2, 0}, {-815.8706, -55.0236, 3}};
    MatrixType result = RealSchur(mat);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (!floatNear(result(i, j), answer(i, j), 1E-5))
                return 1;
    return 0;
}
