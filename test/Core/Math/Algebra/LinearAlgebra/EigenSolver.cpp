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
#include <iostream>
#include "TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"

using namespace Physica::Core;

template<class VectorType>
bool vectorNearZero(const LValueVector<VectorType>& v, double precision) {
    using ScalarType = typename VectorType::ScalarType;
    for (size_t i = 0; i < v.getLength(); ++i)
        if (!scalarNear(v[i], ScalarType::Zero(), precision))
            return false;
    return true;
}

template<class MatrixType>
bool eigenTest(const MatrixType& mat, double precision) {
    EigenSolver solver = EigenSolver(mat, true);
    using ComplexMatrix = typename EigenSolver<MatrixType>::EigenvectorMatrix;
    const size_t order = mat.getRow();
    auto eigenvectors = solver.getEigenvectors();
    for (size_t i = 0; i < order; ++i) {
        auto result = (ComplexMatrix(mat - solver.getEigenvalues()[i] * MatrixType::unitMatrix(order)) * eigenvectors.col(i)).compute();
        if (!vectorNearZero(result[0], precision))
            return false;
    }
    return true;
}

int main() {
    {
        using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3>;
        const MatrixType mat{{-0.590316, -2.19514, -2.37463},
                            {-1.25006, -0.297493, 1.40349},
                            {0.517063, -0.956614, -0.920775}};
        if (!eigenTest(mat, 1E-4))
            return 1;
    }
    return 0;
}
