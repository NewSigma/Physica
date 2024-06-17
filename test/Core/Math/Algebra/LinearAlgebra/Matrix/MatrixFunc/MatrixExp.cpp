/*
 * Copyright 2024 WeiBo He.
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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h>

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseSymmMatrix<ScalarType>;

int main() {
    std::mt19937 gen{};
    const auto mat = MatrixType::random_uniform(16, gen);
    const VectorType v = VectorType::random_uniform(16, gen);

    const auto solver = SymmEigenSolver<MatrixType>(mat, true);
    const VectorType eigenvalues = exp(solver.getEigenvalues());
    const DenseMatrix<ScalarType> eigenvectors = toRealMatrix(solver.getEigenvectors());
    VectorType v1 = eigenvectors.transpose() * v;
    v1 = hadamard(eigenvalues, v1);
    const VectorType answer = eigenvectors * v1;

    const VectorType result = exp(mat.asMatrix()) * v;
    std::cout << answer.format() << '\n' << result.format() << std::endl;
    if (!vectorNear(answer, result, 1E-13))
        return 1;
    return 0;
}
