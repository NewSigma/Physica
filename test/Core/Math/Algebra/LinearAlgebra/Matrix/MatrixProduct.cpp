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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;

int main() {
    using Matrix = DenseMatrix<Scalar<Float, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 4, 4>;
    const Matrix m1{{1, 2}, {2, 1}};
    const Matrix m2{{3, 3}, {1, 5}};
    const Matrix result = m1 * m2;
    const Matrix answer{{9, 9}, {11, 7}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            if (abs((result(i, j) - answer(i, j)) / answer(i, j)).getTrivial() > std::numeric_limits<float>::epsilon())
                return 1;
    return 0;
}
