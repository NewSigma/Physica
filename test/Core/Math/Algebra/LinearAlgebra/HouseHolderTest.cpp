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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"

using namespace Physica::Core;

int main() {
    using T = Scalar<Double, false>;
    {
        const Vector<T> x{2, 3, 4, 5};
        const size_t rank = x.getLength();
        Vector<T> v(rank);
        const T norm = householder(x, v);

        Vector<T> copy = v;
        const T tau = copy[0];
        const T beta = x[0].isNegative() ? norm : -norm;
        copy[0] = 1;

        Vector<T> result = x - tau * (copy * x) * copy;
        if (abs(T(result[0] - beta) / beta) > T(1E-15))
            return 1;
        for (size_t i = 1; i < rank; ++i)
            if (abs(result[i]) > T(1E-14)) //In debug mode, precision can reach 10^-15
                return 1;
    }

    {
        const Vector<T> x{2, 3, 4, 5};
        const size_t rank = x.getLength();
        Vector<T> v(rank);
        const T norm = householder(x, v);

        using MatrixType = DenseMatrix<T, DenseMatrixOption::Column | DenseMatrixOption::Vector, 4, 4>;
        const MatrixType m{x, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        const MatrixType l_answer{{-7.34849, 0, 0, 0}, {-13.0639, 0.203133, -0.729156, -1.66145}, {-20.6846, 0.473976, -1.70137, -3.87671}, {-28.3052, 0.74482, -2.67357, -6.09197}};
        MatrixType l_result = m;

        applyHouseholder(v, l_result);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (!scalarNear(l_result(i, j), l_answer(i, j), 1E-5))
                    return 1;

        const MatrixType r_answer{{-16.3299, -18.2351, -20.1402, -22.0454}, {-0.882225, -0.814514, -0.746803, -0.679092}, {1.15703, 0.913982, 0.67093, 0.427878}, {3.19629, 2.64248, 2.08866, 1.53485}};
        MatrixType r_result = m;
        applyHouseholder(r_result, v);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (!scalarNear(r_result(i, j), r_answer(i, j), 1E-5))
                    return 1;
    }
    /* Complex test */ {
        using ScalarType = ComplexScalar<T>;
        using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector, 2, 2>;

        const Vector<ScalarType> x{{1, 1}, {3, -5}};
        const size_t rank = x.getLength();
        Vector<ScalarType> v(rank);
        const T norm = householder(x, v);

        const MatrixType m{x, {{-2, 7}, {1, 6}}};
        MatrixType householderMat = MatrixType::unitMatrix(2);
        applyHouseholder(v, householderMat);
        MatrixType m1 = householderMat * m;
        MatrixType m2 = m;
        applyHouseholder(v, m2);

        if (!scalarNear(m1(0, 0).norm(), norm, 1E-15))
            return 1;
        if (!scalarNear(m1(1, 0).norm(), T(0), 1E-15))
            return 1;

        m1(1, 0) = m2(1, 0) = T(0);
        if (!matrixNear(m1, m2, 1E-14))
            return 1;
    }
    {
        const Vector<T> x{0, 0, 0, 0};
        Vector<T> v(4);
        const T norm = householder(x, v);
        if (!norm.isZero())
            return 1;
        if (!v.norm().isZero())
            return 1;
    }
    return 0;
}
