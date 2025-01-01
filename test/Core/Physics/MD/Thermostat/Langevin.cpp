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
#include <algorithm>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Physics/MD/Thermostat/Langevin.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

using namespace Physica;
using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, Dynamic>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using ThermoType = Langevin<ScalarType, 1, Dynamic>;
using KineticModel = HardCore<ScalarType, true, Dynamic, RPMDIntegrator::Exact>;
using RandomType = Random<MT19937, 15502868121535481991UL>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double temperatureT = 2;
constexpr double thermostatTime = 0.01;
constexpr double unitMassM = 1;
constexpr size_t numReplica = 8;
constexpr size_t maxHandleNum = 100;
constexpr size_t numSystem = 8;
constexpr size_t numStep = 20000;

MDCellType makeSystem() {
    MDCellType::LatticeMatrix lattice{latticeSize};

    auto posVec = VectorND<ScalarType>::random_uniform<RandomType>(numMolecular);
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
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, numReplica, maxHandleNum);
    ThermoType thermo(temperatureT, thermostatTime, false);

    MDType rpmd = MDType(makeSystem(), numReplica, numReplica, temperatureT, timeStep);
    rpmd.initMomentum<KineticModel, RandomType>();
    kineticModel.updateMass(rpmd.getRingPolymer());

    MatrixType meanCorr(numMolecular, numReplica);
    MatrixType varCorr(numMolecular, numReplica);
    MatrixType temp(numMolecular, numReplica);
    ScalarType meanTemperature = 0;
    ScalarType varTemperature = 0;
    const ScalarType factor = reciprocal(ScalarType(temperatureT * numReplica));
    for (size_t sys = 0; sys < numSystem; ++sys) {
        MatrixType buffer(numMolecular, numReplica, 0);
        ScalarType temperature_sample = 0;
        for (size_t i = 0; i < numStep; ++i) {
            ForceModel forceModel{};
            rpmd.nvt_step<ThermoType, RandomType, KineticModel, ForceModel, SequentialExecutor>(thermo, kineticModel, forceModel);
            auto momentum = rpmd.getRingPolymer().asMatrix().topRows(numMolecular);
            for (size_t replica = 0; replica < numReplica; ++replica) {
                auto col = temp.col(replica);
                col = hadamard(hadamard(momentum.col(replica), momentum.col((replica + 1) % numReplica)), kineticModel.getRepMass()) * factor;
            }
            toNextMean(buffer, i, temp);
            toNextMean(temperature_sample, i, rpmd.calcTemperature<KineticModel>());
        }
        toNextVariance(varCorr, meanCorr, sys, buffer);
        toNextVariance(varTemperature, meanTemperature, sys, temperature_sample);
    }
    const MatrixType deviaCorr = sqrt_elem(varCorr);

    for (size_t i = 0; i < numReplica; ++i)
        for (size_t j = 0; j < numMolecular; ++j)
            if (abs(meanCorr(j, i)) > deviaCorr(j, i) * ScalarType(2.0))
                return 1;
    if (abs(ScalarType(temperatureT) - meanTemperature) > ScalarType(2.0) * sqrt(varTemperature))
        return 1;
    return 0;
}
