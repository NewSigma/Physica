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
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/KSpinRepr.h"

using namespace Physica::Core;
using RealType = float64;
using ScalarType = Complex<RealType>;
using VectorType = VectorND<RealType>;
using RandomType = Random<std::mt19937, 10000>;
constexpr unsigned int NumSite = 6;
constexpr unsigned int NumParticle = NumSite / 2;
constexpr double HoppingT = 1.0;
constexpr double RepelU = 4.0;

void testComplete() {
    const SpinRepr<1, NumSite, true> rRepr(NumParticle, NumParticle);
    const size_t rNumState = rRepr.getNumState();
    size_t kNumState = 0;
    for (unsigned int i = 0; i < NumSite; ++i) {
        const size_t temp = KSpinRepr<1, NumSite, true>(rRepr, i).getNumState();
        kNumState += temp;
    }
    if (rNumState != kNumState)
        exit(EXIT_FAILURE);
}

void testEigen() {
    RealType answer;
    LatticeModel<1> lattice({NumSite}, 1);
    Hubbard<RealType, 1> hubbard(lattice, HoppingT, RepelU);
    {
        using ReprType = SpinRepr<1, NumSite, true>;
        ReprType repr(NumParticle, NumParticle);
        HubbardMatrix<RealType, ReprType> model(hubbard, std::move(repr));

        const size_t numState = model.getNumState();
        JacobiDavidson<RealType> jd(numState, 4);
        jd.compute(model, VectorType::random_uniform(numState, RandomType::getInstance()));
        jd.sort();
        answer = jd.getEigenvalues()[0];
    }
    RealType result;
    {
        using ReprType = KSpinRepr<1, NumSite, true>;
        ReprType repr({NumParticle, NumParticle}, 0);
        HubbardMatrix<ScalarType, ReprType> model(hubbard, std::move(repr));

        const size_t numState = model.getNumState();
        JacobiDavidson<ScalarType> jd(numState, 4);
        jd.compute(model, VectorType::random_uniform(numState, RandomType::getInstance()));
        jd.sort();
        result = jd.getEigenvalues()[0].real();
    }
    if (!scalarNear(answer, result, 1E-14))
        exit(EXIT_FAILURE);
}

int main() {
    testComplete();
    testEigen();
    return 0;
}
