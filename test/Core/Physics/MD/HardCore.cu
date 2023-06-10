/*
 * Copyright 2023 WeiBo He.
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
#include <fstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "Physica/Utils/BenchmarkHelper.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Parallel/StreamPool.cuh"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.cuh"

using namespace Physica::Core;
using namespace Physica::Utils;
constexpr double timeStep = 0.1;
constexpr double collideFactor = 0.001;
const size_t numMolecular = 512;
constexpr double latticeSize = 512;
constexpr double temperatureT = 2;
constexpr double energy = numMolecular * temperatureT / 2;

using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, ScalarType, 1>;

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
        massVec[i] = i % 2 == 0 ? 3.0 : 1.0;
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

void scaleVelocity(MDType& rpmd) {
    const ScalarType energy1 = rpmd.getRingPolymer().calcClassicalKinetic();
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
    using KineticModel = HardCore<ScalarType, 1, CudaExecutor>;
    std::mt19937 gen{};
    MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 100);
    kineticModel.updateMass(rpmd.getRingPolymer());
    rpmd.initMomentum(gen);

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

void testCpuGpuCompare() {
    constexpr int NumData = 11;
    constexpr double precision = 1E-13;
    ScalarType cpu_data[NumData];
    {
        using KineticModel = HardCore<ScalarType, 1>;
        std::mt19937 gen{};
        MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, 100);
        kineticModel.updateMass(rpmd.getRingPolymer());
        rpmd.initMomentum(gen);

        rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, ForceModel());
        cpu_data[0] = calcThermoFlux(rpmd);
        rpmd.nve_step_for<KineticModel, ForceModel, SequentialExecutor>(1.0, kineticModel, ForceModel());
        for (size_t j = 0; j < 10; ++j) {
            cpu_data[j + 1] = calcThermoFlux(rpmd);
            rpmd.nve_step_for<KineticModel, ForceModel, SequentialExecutor>(1.0, kineticModel, ForceModel());
            scaleVelocity(rpmd);
        }
    }
    ScalarType gpu_data[NumData];
    {
        using KineticModel = HardCore<ScalarType, 1, CudaExecutor>;
        std::mt19937 gen{};
        MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 100);
        kineticModel.updateMass(rpmd.getRingPolymer());
        rpmd.initMomentum(gen);

        kineticModel.nve_step(rpmd.getRingPolymer(), timeStep);
        gpu_data[0] = calcThermoFlux(rpmd);
        kineticModel.nve_step_for(1.0, rpmd.getRingPolymer(), timeStep);
        for (size_t j = 0; j < 10; ++j) {
            gpu_data[j + 1] = calcThermoFlux(rpmd);
            kineticModel.do_nve_step_for(1.0, timeStep);
            kineticModel.post_nve_step(rpmd.getRingPolymer());
            scaleVelocity(rpmd);
            kineticModel.updateMomentum(rpmd.getRingPolymer());
        }
    }
    for (int i = 0; i < NumData; ++i) {
        std::cout << gpu_data[i] << ' ' << cpu_data[i] << std::endl;
        if (!scalarNear(gpu_data[i], cpu_data[i], precision))
            exit(EXIT_FAILURE);
    }

}

int main() {
    testMergeStep();
    testCpuGpuCompare();
    return 0;
}
