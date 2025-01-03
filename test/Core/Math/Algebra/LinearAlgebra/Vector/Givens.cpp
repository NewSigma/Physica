/*
 * Copyright 2021-2024 Weibo He.
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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"

using namespace Physica;
using RealType = float64;
using ComplexType = Complex<RealType>;

template<Scalar T>
void test1() {
    Vector2D<T> v{2, 1};
    auto givens_vector = givens(v, 0, 1);
    DenseMatrix<T> v_mat = v;
    applyGivens(givens_vector, v_mat, 0, 1);
    if (abs(v_mat(1, 0).value()) > RealType(1E-15))
        exit(EXIT_FAILURE);
}

int main() {
    test1<RealType>();
    test1<Diff<RealType, DiffMode::Forward, 1>>();
    {
        Vector2D<ComplexType> v{{2, 1}, {1, -3}};
        auto givens_vector = givens(v, 0, 1);
        DenseMatrix<ComplexType> v_mat = v;
        applyGivens(givens_vector, v_mat, 0, 1);
        if (v_mat(1, 0).norm() > RealType(1E-15))
            return 1;
    }
    {
        Vector2D<ComplexType> v{{0.8699464447, 0.1883214037}, {-0.520340944, 0.1297693695}};
        auto givens_vector = givens(v, 0, 1);
        DenseMatrix<ComplexType> v_mat = v;
        applyGivens(givens_vector, v_mat, 0, 1);
        if (v_mat(1, 0).norm() > RealType(1E-15))
            return 1;
    }
    return 0;
}
