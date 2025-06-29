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
#include <iostream>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Parallel/Parallel.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, Dynamic>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, Dynamic, RPMDIntegrator::Exact>;
using RandomSource = Random<MT19937>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double temperatureT = 2;
constexpr double unitMassM = 1;
constexpr size_t numReplica = 8;
constexpr size_t maxHandleNum = 100;

MDCellType makeSystem() {
    MDCellType::LatticeMatrix lattice{latticeSize};

    auto posVec = VectorND<ScalarType>::random_uniform<RandomSource>(numMolecular);
    posVec *= ScalarType(latticeSize);
    std::sort(posVec.begin(), posVec.end());
    MDCellType::PositionMatrix pos(numMolecular, 1);
    pos.col(0) = posVec;

    MDCellType::MassVector massVec(numMolecular);
    for (size_t i = 0; i < numMolecular; ++i) {
        massVec[i] = (i % 2U == 0) ? unitMassM : (unitMassM * 10);
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

int main() {
    const double timeStep = timeStepLambda * (latticeSize / numMolecular) * std::sqrt(unitMassM / temperatureT);
    VectorType nve(20000);
    {
        MDType rpmd = MDType(makeSystem(), numReplica, numReplica, temperatureT, timeStep);
        rpmd.initMomentum<KineticModel, RandomSource>();
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, numReplica, maxHandleNum);
        kineticModel.updateMass(rpmd.getRingPolymer());
        ForceModel forceModel{};
        for (size_t i = 0; i < nve.getLength(); ++i) {
            rpmd.nve_step<Sequential>(kineticModel, forceModel);
            nve[i] = (rpmd.calcClassicalInternalEnergy(forceModel)) / ScalarType(numMolecular * numReplica);
        }
    }
    const bool isEnergyConserved = scalarNear(nve.max(), nve.min(), 1E-11);
    if (!isEnergyConserved)
        return 1;
    return 0;
}
