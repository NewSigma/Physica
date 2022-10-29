/*
 * Copyright 2022 WeiBo He.
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

using namespace Physica::Core;
using ScalarType = Scalar<Float, false>;
using MatrixType = DenseMatrix<ScalarType>;

int main() {
    std::mt19937 gen{};
    const MatrixType A = MatrixType::random_uniform(16, 16, gen);
    const auto d_A = A.toDevice();
    const MatrixType B = d_A.toHost();
    if (A != B)
        return 1;
    return 0;
}
