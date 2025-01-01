/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseSymmMatrix<ScalarType>;
using RandomType = Random<MT19937, std::mt19937::default_seed>;

int main() {
    const auto mat = MatrixType::random_uniform<RandomType>(16);
    const VectorType v = VectorType::random_uniform<RandomType>(16);

    const auto solver = SymmEigenSolver<ScalarType>(mat, true);
    const VectorType eigenvalues = exp(solver.getEigenvalues());
    const DenseMatrix<ScalarType> eigenvectors = solver.getEigenvectors().reals();
    VectorType v1 = eigenvectors.transpose() * v;
    v1 = hadamard(eigenvalues, v1);

    const VectorType answer = eigenvectors * v1;
    const VectorType result = exp(mat) * v;
    if (!vectorNear(answer, result, 1E-13))
        return 1;
    return 0;
}
