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
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Physics/ManyBody/Hubbard1D.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using ReprType = SpinRepr<false>;
constexpr unsigned int NumSite = 3;
constexpr unsigned int NumSpinUp = 1;
constexpr unsigned int NumSpinDown = 1;
constexpr double HoppingT = 1.0;
constexpr double RepelU = 2;

int main() {
    ReprType repr(NumSite, NumSpinUp, NumSpinDown);
    const Hubbard1D<ScalarType, ReprType> model({{NumSite}, 1}, std::move(repr), HoppingT, RepelU);
    const size_t numState = model.getSize();
    MatrixType mat(numState, numState);
    for (size_t i = 0; i < numState; ++i) {
        VectorType temp(numState, 0);
        temp[i] = ScalarType(1);
        auto col = mat.col(i);
        col = model * temp;
    }

    if (!matrixNear(model, mat, 1E-15))
        return 1;
    return 0;
}
