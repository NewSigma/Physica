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
#include "Physica/Core/Physics/ManyBody/Hubbard.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/KSpinRepr.h"

using namespace Physica::Core;
constexpr double HoppingT = 1.0;
constexpr double RepelU = 2;

template<unsigned int NumSite, unsigned int NumSpinUp, unsigned int NumSpinDown>
void testRSpinMatrix1D() {
    using ScalarType = Scalar<Double>;
    using VectorType = Vector<ScalarType>;
    using MatrixType = DenseMatrix<ScalarType>;
    using ReprType = SpinRepr<1, NumSite, NumSpinUp == NumSpinDown>;

    ReprType repr(NumSpinUp, NumSpinDown);
    const Hubbard<ScalarType, ReprType> model({NumSite}, 1, std::move(repr), HoppingT, RepelU);
    const size_t numState = model.getNumState();
    MatrixType mat(numState, numState);
    for (size_t i = 0; i < numState; ++i) {
        VectorType temp(numState, 0);
        temp[i] = ScalarType(1);
        auto col = mat.col(i);
        col = model * temp;
    }

    if (!matrixNear(model, mat, 1E-15))
        exit(EXIT_FAILURE);
}

template<unsigned int NumSiteX, unsigned int NumSiteY, unsigned int NumSpinUp, unsigned int NumSpinDown>
void testRSpinMatrix2D() {
    using ScalarType = Scalar<Double>;
    using VectorType = Vector<ScalarType>;
    using MatrixType = DenseMatrix<ScalarType>;
    using ReprType = SpinRepr<2, NumSiteX * NumSiteY, NumSpinUp == NumSpinDown>;

    ReprType repr(NumSpinUp, NumSpinDown);
    const Hubbard<ScalarType, ReprType> model({NumSiteX, NumSiteY}, 1, std::move(repr), HoppingT, RepelU);
    const size_t numState = model.getNumState();
    MatrixType mat(numState, numState);
    for (size_t i = 0; i < numState; ++i) {
        VectorType temp(numState, 0);
        temp[i] = ScalarType(1);
        auto col = mat.col(i);
        col = model * temp;
    }

    if (!matrixNear(model, mat, 1E-15))
        exit(EXIT_FAILURE);
}

void testKSpinMatrix() {
    constexpr unsigned int NumSite = 4;
    constexpr unsigned int NumParticle = NumSite / 2;
    using ScalarType = ComplexScalar<Scalar<Double>>;
    using VectorType = Vector<ScalarType>;
    using MatrixType = DenseMatrix<ScalarType>;
    using ReprType = KSpinRepr<1, NumSite, true>;

    ReprType repr({NumParticle, NumParticle}, 0);
    Hubbard<ScalarType, ReprType> model({NumSite}, 1, std::move(repr), HoppingT, 4);
    const size_t numState = model.getNumState();
    MatrixType mat(numState, numState);
    for (size_t i = 0; i < numState; ++i) {
        VectorType temp(numState, 0);
        temp[i] = ScalarType(1);
        auto col = mat.col(i);
        col = model * temp;
    }

    if (!matrixNear(model, mat, 1E-15))
        exit(EXIT_FAILURE);
}

int main() {
    testRSpinMatrix1D<3, 2, 1>();
    testRSpinMatrix1D<2, 1, 1>();
    testRSpinMatrix2D<2, 3, 2, 3>();
    testKSpinMatrix();
    return 0;
}
