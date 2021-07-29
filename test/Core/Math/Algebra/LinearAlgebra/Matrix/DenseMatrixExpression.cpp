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
    DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3> mat1{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
    DenseMatrix<Scalar<Float, false>, DenseMatrixOption::Column | DenseMatrixOption::Element, 3, 3> mat2{1, 1, 1, 1, 1, 1, 1, 1, 1};
    {
        DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Row | DenseMatrixOption::Vector, 3, 3> mat = -(mat1 + mat2);
        for (auto ite = mat.cbegin(); ite != mat.cend(); ++ite)
            for (auto ite1 = (*ite).cbegin(); ite1 != (*ite).cend(); ++ite1)
                if (*ite1 != Scalar<Double, false>(-2))
                    return 1;
    }
    {
        DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Row | DenseMatrixOption::Vector, 3, 3> mat = mat1 * mat2;
        for (auto ite = mat.cbegin(); ite != mat.cend(); ++ite)
            for (auto ite1 = (*ite).cbegin(); ite1 != (*ite).cend(); ++ite1)
                if (*ite1 != Scalar<Double, false>(3))
                    return 1;
    }
    return 0;
}
