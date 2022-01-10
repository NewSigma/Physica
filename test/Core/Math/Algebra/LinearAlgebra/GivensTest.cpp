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
#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;

int main() {
    using RealType = Scalar<Double, false>;
    using ComplexType = ComplexScalar<RealType>;
    {
        Vector<RealType, 2> v{2, 1};
        auto givens_vector = givens(v, 0, 1);
        DenseMatrix<RealType> v_mat = v;
        applyGivens(givens_vector, v_mat, 0, 1);
        if (abs(v_mat(1, 0)) > RealType(1E-15))
            return 1;
    }
    {
        Vector<ComplexType, 2> v{{2, 1}, {1, -3}};
        auto givens_vector = givens(v, 0, 1);
        DenseMatrix<ComplexType> v_mat = v;
        applyGivens(givens_vector, v_mat, 0, 1);
        if (v_mat(1, 0).norm() > RealType(1E-15))
            return 1;
    }
    return 0;
}
