/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseSymmMatrix<ScalarType>;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

int main() {
    const auto mat = MatrixType::random_uniform<RandomSource>(16);
    const VectorType v = VectorType::random_uniform<RandomSource>(16);

    const auto solver = SymmEigenSolver<ScalarType>(mat, true);
    const VectorType eigenvalues = exp(solver.getEigenvalues());
    const DenseMatrix<ScalarType> eigenvectors = solver.getEigenvectors().reals();
    VectorType v1 = eigenvectors.transpose() * v;
    v1 = hadamard(eigenvalues, v1);

    const VectorType answer = eigenvectors * v1;
    const VectorType result = exp(mat) * v;
    expect(vectorNear(answer, result, 1E-13));
    return 0;
}
