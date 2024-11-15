/*
 * Copyright 2021 Weibo He.
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
#include <iostream>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>
#include <Physica/Core/MultiPrecision/Diff.h>

using namespace Physica::Core;

int main() {
    {
        using Matrix = DenseMatrix<float32, MatrixOption::Col | MatrixOption::Vector, 2, 2>;
        const Matrix m1{{1, 2}, {2, 1}};
        const Matrix m2{{3, 3}, {1, 5}};
        const Matrix result = m1 * m2;
        const Matrix answer{{9, 9}, {11, 7}};
        if (!matrixNear(result, answer, std::numeric_limits<float>::epsilon()))
            return 1;
    }
    {
        using ScalarType = Diff<float32, DiffMode::Reverse>;
        using MatrixType = DenseMatrix<float32, MatrixOption::Col | MatrixOption::Element, 3, 3>;
        using DVector = Vector3D<ScalarType>;
        using DMatrix = DenseMatrix<ScalarType, MatrixOption::Col | MatrixOption::Element, 3, 3>;
        const DMatrix m{1, 2, 3, 4, 5, 6, 7, 8, 9};
        const DVector x{1, 2, 3};
        const DVector y = m * x;
        y.sum().reverse();

        const Vector3D<float32> dx = toGradVector(x);
        for (size_t i = 0; i < dx.getLength(); ++i)
            if (!scalarNear(dx[i], toValueMatrix(m).col(i).sum(), 1E-15))
                return 1;

        const MatrixType dm = toGradMatrix(m);
        for (size_t r = 0; r < dm.getRow(); ++r)
            for (size_t c = 0; c < dm.getCol(); ++c)
                if (!scalarNear(dm(r, c), toValueVector(x).calc(c), 1E-15))
                    return 1;
    }
    return 0;
}
