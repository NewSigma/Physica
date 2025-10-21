/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

namespace {
    void test1() {
        using MatrixType = DenseMatrix<cfloat64, MatrixOption::Row>;
        const MatrixType data = MatrixType::random_uniform<RandomSource>(64);
        const DenseHermiteMatrix<cfloat64> hermite = data + data.hermite();

        EigenSolver<cfloat64> eig(hermite, false);
        eig.sort();
        JacobiDavidson<cfloat64> jd(hermite.getRow(), 48);
        jd.compute(hermite, VectorND<cfloat64>::random_uniform<RandomSource>(data.getRow()));
        jd.sort();

        if (!vectorNear(jd.getEigenvalues(), eig.getEigenvalues().head(jd.getNumRequired()), 1E-13))
            exit(1);
    }
}

int main() {
    test1();
    return 0;
}
