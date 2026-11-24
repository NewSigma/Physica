/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Parallel/Parallel.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, Dynamic>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, Dynamic, RPMDIntegrator::Exact>;
using RandomSource = Random<>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double temperatureT = 2;
constexpr double unitMassM = 1;
constexpr size_t numReplica = 8;
constexpr size_t maxHandleNum = 100;

namespace {
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
    expect<RandomSource>(isEnergyConserved);
    return 0;
}
