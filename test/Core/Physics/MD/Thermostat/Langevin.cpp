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
#include <algorithm>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Physics/MD/Thermostat/Langevin.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Parallel/Parallel.h"
#include "Test.h"

using namespace Physica;
using namespace Physica;
using T = float64;
using MDType = RPMD<T, 1, Dynamic>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<T, 1>;
using ThermoType = Langevin<T, 1, Dynamic>;
using KineticModel = HardCore<T, true, Dynamic, RPMDIntegrator::Exact>;
using RandomSource = Random<MT19937, 15502868121535481991UL>;
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

namespace {
    MDCellType makeSystem() {
        MDCellType::LatticeMatrix lattice{latticeSize};

        auto posVec = VectorND<T>::random_uniform<RandomSource>(numMolecular);
        std::ranges::sort(posVec.begin(), posVec.end());
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
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, numReplica, maxHandleNum);
    ThermoType thermo(temperatureT, thermostatTime, false);

    MDType rpmd = MDType(makeSystem(), numReplica, numReplica, temperatureT, timeStep);
    rpmd.initMomentum<KineticModel, RandomSource>();
    kineticModel.updateMass(rpmd.getRingPolymer());

    MatrixND<T> meanCorr(numMolecular, numReplica);
    MatrixND<T> varCorr(numMolecular, numReplica);
    MatrixND<T> temp(numMolecular, numReplica);
    T meanTemperature = 0;
    T varTemperature = 0;
    const T factor = reciprocal(T(temperatureT * numReplica));
    for (size_t sys = 0; sys < numSystem; ++sys) {
        MatrixND<T> buffer(numMolecular, numReplica, 0);
        T temperature_sample = 0;
        for (size_t i = 0; i < numStep; ++i) {
            ForceModel forceModel{};
            rpmd.nvt_step<RandomSource, Sequential>(thermo, kineticModel, forceModel);
            auto momentum = rpmd.getRingPolymer().asMatrix().topRows(numMolecular);
            for (size_t replica = 0; replica < numReplica; ++replica) {
                auto col = temp.col(replica);
                col = hadamard(hadamard(momentum.col(replica), momentum.col((replica + 1) % numReplica)), kineticModel.getRepMass()) * factor;
            }
            buffer.toNextMean(i, temp);
            temperature_sample.toNextMean(i, rpmd.calcTemperature<KineticModel>());
        }
        varCorr.toNextVariance(meanCorr, sys, buffer);
        varTemperature.toNextVariance(meanTemperature, sys, temperature_sample);
    }
    const MatrixND<T> deviaCorr = sqrt_elem(varCorr);

    for (size_t i = 0; i < numReplica; ++i)
        for (size_t j = 0; j < numMolecular; ++j)
            expect(abs(meanCorr[j, i]) < deviaCorr[j, i] * T(2.0));

    expect(abs(T(temperatureT) - meanTemperature) < T(2.0) * sqrt(varTemperature));
    return 0;
}
