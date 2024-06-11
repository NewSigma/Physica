/*
 * Copyright 2022-2024 WeiBo He.
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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h>

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using ComplexType = ComplexScalar<ScalarType>;

int main() {
    using MatrixType = DenseMatrix<ComplexType, MatrixOption::Row | MatrixOption::Vector>;
    std::mt19937 gen{};
    const MatrixType data = MatrixType::random_uniform(64, gen);
    const DenseHermiteMatrix<ComplexType> hermite = data + data.hermite();

    EigenSolver<MatrixType> eig(hermite, false);
    eig.sort();
    JacobiDavidson<ComplexType> jd(hermite.getRow(), 61);
    jd.compute(hermite, Vector<ComplexType>::random_uniform(data.getRow(), gen));
    jd.sort();

    if (!vectorNear(jd.getEigenvalues(), eig.getEigenvalues().head(jd.getNumRequired()), 1E-13))
        return 1;
    return 0;
}
