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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica::Core;
using ScalarType = Real<Double>;
using ComplexType = Complex<ScalarType>;
using MatrixType = DenseHermiteMatrix<ComplexType>;

int main() {
    auto& gen = Random<MT19937, std::mt19937::default_seed>::getInstance();
    const auto mat = MatrixType::random_uniform(5000, gen);

    JacobiDavidson<ComplexType> jd(mat.getRow(), 4);
    jd.compute(mat, VectorND<ComplexType>::random_uniform(mat.getRow(), gen));
    jd.sort();

    std::cout << toRealVector(jd.getEigenvalues()).format() << std::endl;
    return 0;
}
