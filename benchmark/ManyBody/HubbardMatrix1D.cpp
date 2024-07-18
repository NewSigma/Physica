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
#include <gperftools/profiler.h>
#include <Physica/Core/Math/Random/RandomPool.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h>
#include <Physica/Core/Physics/ManyBody/ExactDiag/Hamilton/HubbardMatrix.h>
#include <Physica/Core/Physics/ManyBody/ExactDiag/ReprSpace/KSpinRepr.h>
#include <Physica/Core/Parallel/Executor/ThreadExecutor.h>
#include <Physica/Utils/Cycler.h>

using namespace Physica::Core;
using namespace Physica::Utils;
using RealType = Scalar<Double>;
using ScalarType = ComplexScalar<RealType>;
using VectorType = Vector<RealType>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr unsigned int NumSite = 10;
constexpr unsigned int NumParticle = NumSite / 2;
constexpr double HoppingT = 1.0;
constexpr double RepelU = 4;

int main() {
    const size_t kSize = FFT<RealType, 1>::rSizeToKSize(NumSite);
    //ProfilerStart("profiler.dat");
    const LatticeModel<1> lattice({NumSite}, 1);
    const Hubbard<RealType, 1> hubbard(lattice, HoppingT, RepelU);
    for (unsigned int kIndex = 0; kIndex < kSize; ++kIndex) {
        using ReprType = KSpinRepr<1, NumSite, true>;
        ReprType repr({NumParticle, NumParticle}, kIndex);
        const HubbardMatrix<ScalarType, ReprType> model(hubbard, std::move(repr));

        const size_t numState = model.getNumState();
        JacobiDavidson<ScalarType> jd(numState, 4);
        jd.compute(model, VectorType::random_uniform(numState, RandomPoolType::getInstance().getGen()));
    }
    //ProfilerStop();
    return 0;
}
