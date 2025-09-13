/*
 * Copyright 2023-2025 Weibo He.
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
#include <benchmark/benchmark.h>
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Statistics/ProbDistribution.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using RandomSource = Random<MT19937>;
using PDFType = ProbDistribution<ScalarType>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 512;
constexpr size_t numMolecular = 512;
constexpr double temperatureT = 2;
constexpr double unitMassM = 1;
constexpr size_t numReplica = 1;
constexpr size_t maxHandleNum = 1000;
using MDType = RPMD<ScalarType, 1, numReplica>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, numReplica, RPMDIntegrator::Exact>;

namespace {
    MDCellType makeSystem() {
        MDCellType::LatticeMatrix lattice{latticeSize};

        auto posVec = VectorND<ScalarType>::random_uniform<RandomSource>(numMolecular);
        std::sort(posVec.begin(), posVec.end());
        MDCellType::PositionMatrix pos(numMolecular, 1);
        pos.col(0) = posVec;

        MDCellType::MassVector massVec(numMolecular);
        for (size_t i = 0; i < numMolecular; ++i)
            massVec[i] = unitMassM;
        return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
    }

    void func(benchmark::State& state) {
        const double timeStep = timeStepLambda * (latticeSize / numMolecular) * std::sqrt(unitMassM / temperatureT);
        MDType rpmd = MDType(makeSystem(), numReplica, numReplica, temperatureT, timeStep);
        rpmd.initMomentum<KineticModel, RandomSource>();
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, numReplica, maxHandleNum);
        ForceModel forceModel{};
        kineticModel.updateMass(rpmd.getRingPolymer());

        for (auto _ : state)
            rpmd.nve_step<Sequential>(kineticModel, forceModel);
    }
}

BENCHMARK(func)->Name("HardCore")->Unit(benchmark::kMicrosecond);
