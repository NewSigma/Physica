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
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.cuh"
#include "Physica/Core/Math/Statistics/ProbabilityDistributionFunction.h"

using namespace Physica::Core;
constexpr double timeStep = 0.001;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 50;
constexpr size_t numMolecular = 50;
constexpr double energy = 50;
constexpr size_t numStep = 100000000;

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
        massVec[i] = i % 2 == 0 ? M_SQRT2 : 1.0;
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

void scaleVelocity(MDType& rpmd) {
    const ScalarType energy1 = rpmd.getRingPolymer().calcClassicalKinetic();
    auto phase = rpmd.getPhaseMatrix().col(0);
    auto momentum = phase.head(numMolecular);
    momentum *= sqrt(ScalarType(energy) / energy1);
}

void run(unsigned int sys, std::mt19937& gen, const ProbabilityDistributionFunction<ScalarType>& originPdf, MatrixType& record) {
    MDType rpmd = MDType(makeSystem(gen), 1, 1, 1, timeStep);
    KineticModel kineticModel(latticeSize, collideFactor, numMolecular);
    kineticModel.updateMass(rpmd.getRingPolymer());

    rpmd.initMomentum(gen);
    scaleVelocity(rpmd);

    ProbabilityDistributionFunction<ScalarType> pdf(originPdf.getFromPoint(), originPdf.getToPoint(), originPdf.getNumBin());
    for (size_t i = 0; i < numStep; ++i) {
        rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, ForceModel());
        pdf.sample(rpmd.getPhaseMatrix()(numMolecular / 2, 0));
        if (i % 10000 == 0)
            scaleVelocity(rpmd);
    }
    record[sys] = pdf.makeDistribution();
}
