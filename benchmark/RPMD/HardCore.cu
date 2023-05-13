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
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 32;
constexpr double energy = 32;

using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = FreeModel<ScalarType, ScalarType, 1>;
using KineticModel = HardCore<ScalarType, CudaExecutor>;

MDCellType makeSystem(size_t numMolecular, std::mt19937& gen) {
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
        massVec[i] = i % 2 == 0 ? M_SQRT2 : 1.0;
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

void scaleVelocity(size_t numMolecular, MDType& rpmd) {
    const ScalarType energy1 = rpmd.getRingPolymer().calcClassicalKinetic();
    auto phase = rpmd.getPhaseMatrix().col(0);
    auto momentum = phase.head(numMolecular);
    momentum *= sqrt(ScalarType(energy) / energy1);
}

void run(double timeStep, std::mt19937& gen) {
    const size_t numStep = 10000;
    const size_t numMolecular = 1024;

    MDType rpmd = MDType(makeSystem(numMolecular, gen), 1, 1, 1, timeStep);
    KineticModel kineticModel(latticeSize, collideFactor, numMolecular);
    kineticModel.updateMass(rpmd.getRingPolymer());
    rpmd.initMomentum(gen);
    scaleVelocity(numMolecular, rpmd);

    for (size_t i = 0; i < numStep; ++i) {
        rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, ForceModel());
    }
}
