/*
 * Copyright 2023-2024 Weibo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Core/Math/Random/RandomPool.h"

using namespace Physica;
using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, Dynamic>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, Dynamic, RPMDIntegrator::Exact>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double temperatureT = 2;
constexpr double unitMassM = 1;
constexpr size_t numReplica = 8;
constexpr size_t maxHandleNum = 100;

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
        massVec[i] = (i % 2U == 0) ? unitMassM : (unitMassM * 10);
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

int main() {
    const double timeStep = timeStepLambda * (latticeSize / numMolecular) * std::sqrt(unitMassM / temperatureT);
    auto& pool = RandomPoolType::getInstance();
    auto& gen = pool.getGen();
    VectorType nve(20000);
    {
        MDType rpmd = MDType(makeSystem(gen), numReplica, numReplica, temperatureT, timeStep);
        rpmd.initMomentum<KineticModel, decltype(gen)>(gen);
        KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, numReplica, maxHandleNum);
        kineticModel.updateMass(rpmd.getRingPolymer());
        ForceModel forceModel{};
        for (size_t i = 0; i < nve.getLength(); ++i) {
            rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, forceModel);
            nve[i] = (rpmd.calcClassicalInternalEnergy(forceModel)) / ScalarType(numMolecular * numReplica);
        }
    }
    const bool isEnergyConserved = scalarNear(nve.max(), nve.min(), 1E-11);
    if (!isEnergyConserved)
        return 1;
    return 0;
}
