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
#include <random>
#include <algorithm>
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.cuh"

using namespace Physica::Core;
constexpr double timeStep = 0.1;
constexpr double collideFactor = 0.001;
const size_t numMolecular = 512;
constexpr double latticeSize = 512;
constexpr double temperatureT = 2;
constexpr double energy = numMolecular * temperatureT / 2;

using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using RandomType = Random<MT19937, 12345>;
using MDType = RPMD<ScalarType, 1, 1>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;

MDCellType makeSystem() {
    MDCellType::LatticeMatrix lattice{latticeSize};

    std::uniform_real_distribution dist{};
    VectorND<ScalarType> posVec(numMolecular);
    for (auto& elem : posVec)
        elem = dist(RandomType::getInstance()) * latticeSize;
    std::sort(posVec.begin(), posVec.end());
    MDCellType::PositionMatrix pos(numMolecular, 1);
    pos.col(0) = posVec;

    MDCellType::MassVector massVec(numMolecular);
    for (size_t i = 0; i < numMolecular; ++i) {
        massVec[i] = i % 2 == 0 ? 3.0 : 1.0;
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

void scaleVelocity(MDType& rpmd) {
    const ScalarType energy1 = rpmd.getRingPolymer().calcKineticClassical();
    auto phase = rpmd.getPhaseMatrix().col(0);
    auto momentum = phase.head(numMolecular);
    momentum *= sqrt(ScalarType(energy) / energy1);
}

ScalarType calcThermoFlux(MDType& rpmd) {
    ScalarType flux = 0;
    auto col = rpmd.getPhaseMatrix().col(0);
    const auto& massVec = rpmd.getMassVec();
    for (size_t i = 0; i < numMolecular; ++i)
        flux += col[i] * square(col[i]) / square(massVec[i]);
    return flux;
}

void testMergeStep() {
    using KineticModel = HardCore<ScalarType, true, 1, RPMDIntegrator::Exact, CUDAExecutor>;
    MDType rpmd = MDType(makeSystem(), 1, 1, temperatureT, timeStep);
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 100);
    kineticModel.updateMass(rpmd.getRingPolymer());
    rpmd.initMomentum<KineticModel, RandomType>();

    const MatrixType origin = rpmd.getPhaseMatrix();
    kineticModel.nve_step(rpmd.getRingPolymer(), timeStep);
    kineticModel.nve_step(rpmd.getRingPolymer(), timeStep);
    const MatrixType mat1 = rpmd.getPhaseMatrix();

    rpmd.getPhaseMatrix() = origin;
    kineticModel.pre_nve_step(rpmd.getRingPolymer());
    kineticModel.do_nve_step(timeStep, 2);
    kineticModel.post_nve_step(rpmd.getRingPolymer());
    const MatrixType mat2 = rpmd.getPhaseMatrix();
    if (!matrixNear(mat1, mat2, 1E-5)) {
        printf("testSeperateStep failed\n");
        exit(EXIT_FAILURE);
    }
}

template<bool IsFixedBoundary>
void testCpuGpuCompare() {
    constexpr int NumData = 8;
    constexpr double precision = 1E-13;
    ScalarType cpu_data[NumData];
    {
        using KineticModel = HardCore<ScalarType, IsFixedBoundary, 1, RPMDIntegrator::Exact>;
        std::ignore = RandomType::getInstance().reseed();
        MDType rpmd = MDType(makeSystem(), 1, 1, temperatureT, timeStep);
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, 100);
        kineticModel.updateMass(rpmd.getRingPolymer());
        rpmd.initMomentum<KineticModel, RandomType>();

        ForceModel forceModel{};
        rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, forceModel);
        cpu_data[0] = calcThermoFlux(rpmd);
        rpmd.nve_step_for<KineticModel, ForceModel, SequentialExecutor>(1.0, kineticModel, forceModel);
        for (size_t j = 1; j < NumData; ++j) {
            cpu_data[j] = calcThermoFlux(rpmd);
            rpmd.nve_step_for<KineticModel, ForceModel, SequentialExecutor>(1.0, kineticModel, forceModel);
            scaleVelocity(rpmd);
        }
    }
    ScalarType gpu_data[NumData];
    {
        using KineticModel = HardCore<ScalarType, IsFixedBoundary, 1, RPMDIntegrator::Exact, CUDAExecutor>;
        std::ignore = RandomType::getInstance().reseed();
        MDType rpmd = MDType(makeSystem(), 1, 1, temperatureT, timeStep);
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 100);
        kineticModel.updateMass(rpmd.getRingPolymer());
        rpmd.initMomentum<KineticModel, RandomType>();

        kineticModel.nve_step(rpmd.getRingPolymer(), timeStep);
        gpu_data[0] = calcThermoFlux(rpmd);
        kineticModel.nve_step_for(1.0, rpmd.getRingPolymer(), timeStep);
        for (size_t j = 1; j < NumData; ++j) {
            gpu_data[j] = calcThermoFlux(rpmd);
            kineticModel.do_nve_step_for(1.0, timeStep);
            kineticModel.post_nve_step(rpmd.getRingPolymer());
            scaleVelocity(rpmd);
            kineticModel.updateMomentum(rpmd.getRingPolymer());
        }
    }
    for (int i = 0; i < NumData; ++i) {
        if (!scalarNear(gpu_data[i], cpu_data[i], precision))
            exit(EXIT_FAILURE);
    }
}

int main() {
    testMergeStep();
    testCpuGpuCompare<false>();
    testCpuGpuCompare<true>();
    return 0;
}
