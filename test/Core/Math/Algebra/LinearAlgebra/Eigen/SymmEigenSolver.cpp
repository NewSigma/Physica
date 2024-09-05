/*
 * Copyright 2022-2023 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"

using namespace Physica::Core;

template<class MatrixType>
bool eigenTest(const MatrixType& mat, double precision) {
    using ScalarType = typename MatrixType::ScalarType;
    using RealType = typename ScalarType::RealType;
    using ComplexVector = Vector<ComplexScalar<RealType>, MatrixType::RowAtCompile>;
    using EigenvectorMatrix = typename SymmEigenSolver<ScalarType>::EigenvectorMatrix;

    auto solver = SymmEigenSolver<ScalarType>(mat, true);
    solver.sort();

    const size_t order = mat.getRow();
    EigenvectorMatrix eigenvectors = solver.getEigenvectors();
    for (size_t i = 0; i < order; ++i) {
        if (i > 1 && solver.getEigenvalues()[i - 1] > solver.getEigenvalues()[i])
            return false;
        ComplexVector v1 = mat * eigenvectors.col(i);
        ComplexVector v2 = solver.getEigenvalues()[i] * eigenvectors.col(i).asVector();
        if (!vectorNear(v1, v2, precision))
            return false;
    }
    return true;
}

int main() {
    using RealType = Scalar<Double>;
    {
        using MatrixType = DenseMatrix<RealType, MatrixOption::Column | MatrixOption::Vector, 3, 3>;
        const MatrixType mat1{{-0.590316, -2.19514, -2.37463},
                             {-1.25006, -0.297493, 1.40349},
                             {0.517063, -0.956614, -0.920775}};
        MatrixType mat = mat1 + mat1.transpose();
        if (!eigenTest(mat, 1E-14))
            return 1;
    }
    {
        using MatrixType = DenseSymmMatrix<RealType>;
        std::mt19937 gen{};
        const auto mat = MatrixType::random_uniform(8, gen);
        if (!eigenTest(mat, 1E-14))
            return 1;
    }
}
