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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Test.h"

using namespace Physica;
using RealType = float64;
using ComplexType = Complex<RealType>;

namespace {
    template<Scalar T>
    void test1() {
        Vector2D<T> v{T(2), T(1)};
        auto givens_vector = givens(v, 0, 1);
        DenseMatrix<T> v_mat = v;
        applyGivens(givens_vector, v_mat, 0, 1);
        expect(abs(v_mat[1, 0].value()) < RealType(1E-15));
    }

    void emptyTest() {
        Vector2D<float64> v{1, 0};
        auto g = givens(v, 0, 1);
        expect(g[0] == float64(1));
        expect(g[1] == float64(0));
    }
}

int main() {
    test1<RealType>();
    test1<Diff<RealType, DiffMode::Forward, 1>>();
    emptyTest();
    {
        Vector2D<ComplexType> v{{2, 1}, {1, -3}};
        auto givens_vector = givens(v, 0, 1);
        DenseMatrix<ComplexType> v_mat = v;
        applyGivens(givens_vector, v_mat, 0, 1);
        expect(v_mat[1, 0].norm() < RealType(1E-15));
    }
    {
        Vector2D<ComplexType> v{{0.8699464447, 0.1883214037}, {-0.520340944, 0.1297693695}};
        auto givens_vector = givens(v, 0, 1);
        DenseMatrix<ComplexType> v_mat = v;
        applyGivens(givens_vector, v_mat, 0, 1);
        expect(v_mat[1, 0].norm() < RealType(1E-15));
    }
    return 0;
}
