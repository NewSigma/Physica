/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Physics/MD/Thermostat/Langevin.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Statistics/ProbabilityDistributionFunction.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Utils/BenchmarkHelper.h"

using namespace Physica::Core;
using ScalarType = float64;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using PDFType = ProbabilityDistributionFunction<ScalarType>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 512;
constexpr size_t numMolecular = 512;
constexpr double temperatureT = 2;
constexpr double unitMassM = 1;
constexpr size_t numReplica = 1;
constexpr size_t maxHandleNum = 100;
using MDType = RPMD<ScalarType, 1, numReplica>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, numReplica, RPMDIntegrator::Exact>;

namespace {
    MDCellType makeSystem(std::mt19937& gen) {
        typename MDCellType::LatticeMatrix lattice{latticeSize};

        std::uniform_real_distribution dist{};
        Vector<ScalarType> posVec(numMolecular);
        for (auto& elem : posVec)
            elem = dist(gen) * latticeSize;
        std::sort(posVec.begin(), posVec.end());
        typename MDCellType::PositionMatrix pos(numMolecular, 1);
        pos.col(0) = posVec;

        typename MDCellType::MassVector massVec(numMolecular);
        for (size_t i = 0; i < numMolecular; ++i) {
            massVec[i] = unitMassM;
        }
        return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
    }

    static void main(benchmark::State& state) {
        const double timeStep = timeStepLambda * (latticeSize / numMolecular) * std::sqrt(unitMassM / temperatureT);
        std::mt19937 gen{};
        MDType rpmd = MDType(makeSystem(gen), numReplica, numReplica, temperatureT, timeStep);
        rpmd.initMomentum<KineticModel, decltype(gen)>(gen);
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, numReplica, maxHandleNum);
        kineticModel.updateMass(rpmd.getRingPolymer());

        for (auto _ : state) {
            ForceModel forceModel{};
            for (size_t i = 0; i < 50000; ++i)
                rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, forceModel);
        };
    }
}

BENCHMARK(main)->Name("HardCore")->Unit(benchmark::kSecond);
