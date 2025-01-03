/*
 * Copyright 2022-2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Diff.h"

using namespace Physica;
using RandomType = Random<MT19937, std::mt19937::default_seed>;

template<Matrix T>
bool eigenTest(const T& mat, double precision) {
    using ScalarType = T::ScalarType;
    using VectorType = DenseVector<ScalarType, T::RowAtCompile>;
    using EigenvectorMatrix = SymmEigenSolver<ScalarType>::EigenvectorMatrix;

    auto solver = SymmEigenSolver<ScalarType>(mat, true);
    solver.sort();

    const size_t order = mat.getRow();
    EigenvectorMatrix eigenvectors = solver.getEigenvectors();
    for (size_t i = 0; i < order; ++i) {
        if (i > 1 && solver.getEigenvalues()[i - 1] > solver.getEigenvalues()[i])
            return false;
        VectorType v1 = mat * eigenvectors.col(i);
        VectorType v2 = eigenvectors.col(i) * ScalarType(solver.getEigenvalues()[i]);
        if (!vectorNear(v1, v2, precision))
            return false;
    }
    return true;
}

int main() {
    {
        using MatrixType = DenseMatrix<float64, MatrixOption::Col | MatrixOption::Vector, 3, 3>;
        const MatrixType data{
                {-0.590316, -2.195140, -2.374630},
                {-1.250060, -0.297493,  1.403490},
                { 0.517063, -0.956614, -0.920775}
        };
        const MatrixType mat = data + data.transpose();
        if (!eigenTest(mat, 1E-14))
            return 1;

        using T = Diff<float64, DiffMode::Forward, 1>;
        using DiffMatrix = DenseMatrix<T, MatrixOption::Col | MatrixOption::Vector, 3, 3>;
        const DiffMatrix mat1 = mat;
        if (!eigenTest(mat1, 1E-14))
            return 1;
    }
    {
        using MatrixType = DenseSymmMatrix<float64>;
        const auto mat = MatrixType::random_uniform<RandomType>(8);
        if (!eigenTest(mat, 1E-13))
            return 1;
    }
}
