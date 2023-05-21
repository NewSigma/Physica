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
#include <algorithm>
#include <fstream>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.cuh"

using namespace Physica::Core;
constexpr double timeStep = 0.1;
constexpr double collideFactor = 0.005;
const size_t numMolecular = 512;
constexpr double latticeSize = 512;
constexpr double temperatureT = 2;
constexpr double energy = numMolecular * temperatureT / 2;
constexpr size_t numSample = 10000;

using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = FreeModel<ScalarType, ScalarType, 1>;
using KineticModel = HardCore<ScalarType, CudaExecutor>;

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

void run(unsigned int sys, MatrixType& record, std::mt19937& gen) {
    MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
    KineticModel kineticModel(latticeSize, collideFactor, numMolecular, 100);
    kineticModel.updateMass(rpmd.getRingPolymer());
    rpmd.initMomentum(gen);
    scaleVelocity(rpmd);

    Vector<ScalarType> mean(record.getRow(), 0);
    for (size_t sample = 0; sample < numSample; ++sample) {
        kineticModel.pre_nve_step(rpmd.getRingPolymer());
        kineticModel.do_nve_step_for(1.0, rpmd.getRingPolymer(), timeStep);
        const ScalarType flux0 = calcThermoFlux(rpmd);
        for (size_t j = 0; j < mean.getLength(); ++j) {
            const ScalarType flux = calcThermoFlux(rpmd);
            toNextMean(mean[j], sample, flux0 * flux);
            kineticModel.do_nve_step_for(1.0, rpmd.getRingPolymer(), timeStep);
            scaleVelocity(rpmd);
            kineticModel.pre_nve_step(rpmd.getRingPolymer());
        }
    }
    record[sys] = std::move(mean);
}
