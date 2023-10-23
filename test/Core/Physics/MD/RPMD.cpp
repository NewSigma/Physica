/*
 * Copyright 2022-2023 WeiBo He.
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
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica::Core;
using RandomGenerator = std::mt19937;
using RandomPoolType = RandomPool<RandomGenerator, 3438603950906262893>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

template<class ScalarType>
RPMD<ScalarType, ScalarType> makeSystem(RandomGenerator& gen) {
    using MDCellType = typename RPMD<ScalarType, ScalarType>::MDCellType;
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType, ScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}

bool testDriftMomentum(double precision) {
    using ScalarType = Scalar<Double>;
    auto& gen = RandomPoolType::getGen();
    auto rpmd = makeSystem<ScalarType>(gen);
    rpmd.initMomentum(gen);
    for (int i = 0; i < 3; ++i) {
        ScalarType sum = 0;
        for (size_t j = i; j < rpmd.getDOF(); j += 3)
            sum += rpmd.getPhaseMatrix().row(j).asVector().sum();
        if (!scalarNear(sum, ScalarType(0), precision))
            return false;
    }
    return true;
}

bool testCalcKinetic(double precision) {
    using ScalarType = Scalar<Double>;
    using PosScalarType = ScalarType;
    using ForceModel = SilveraGoldman<ScalarType, PosScalarType>;

    auto& gen = RandomPoolType::getGen();
    auto rpmd = makeSystem<ScalarType>(gen);
    rpmd.initMomentum(gen);
    ForceModel forceModel(pair_cutoff);
    rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);

    const ScalarType kinetic1 = rpmd.calcKinetic();
    ScalarType kinetic2 = 0;
    for (size_t i = 0; i < rpmd.getDOF(); ++i)
        kinetic2 += rpmd.calcKinetic(i);
    return scalarNear(kinetic1, kinetic2, precision);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
void testMDRun() {
    using ScalarType = Scalar<Double>;
    using PosScalarType = ScalarType;
    using ForceModel = SilveraGoldman<ScalarType, PosScalarType>;
    using ThermostatType = DoubleThermo<ScalarType, PosScalarType>;
    using KineticModel = FreeModel<ScalarType, PosScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
    ScalarType mean = 0;
    ScalarType var = 0;
    {
        const ThermostatType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        ForceModel forceModel(pair_cutoff);
        RandomGenerator& gen = RandomPoolType::getGen();
        auto rpmd = makeSystem<ScalarType>(gen);
        rpmd.initMomentum(gen);
        rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);

        for (unsigned int i = 0; i < 6; ++i) {
            ScalarType temp = 0;
            rpmd.nvt_step_for<ThermostatType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-12),
                thermo,
                kineticModel,
                forceModel);

            for (unsigned int j = 0; j < 100; ++j) {
                rpmd.nvt_step<ThermostatType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(thermo, kineticModel, forceModel);
                toNextMean(temp, j, rpmd.calcKinetic());
            }
            toNextVariance(var, mean, i, temp);
        }
    }
    constexpr double answer = 61.8;
    const ScalarType energyPerMol = PhyConst<AU>::temperatureToK(double(mean) / numMolecular);
    const ScalarType deviation = PhyConst<AU>::temperatureToK(std::sqrt(double(var))) / numMolecular;
    if (abs(energyPerMol - answer) > deviation * 2.0)
        exit(EXIT_FAILURE);
}

int main() {
    ThreadPool::numThreadRequired = 4;
    testDriftMomentum(1E-12);
    testCalcKinetic(1E-14);
    testMDRun();
    ThreadPool::getInstance().shouldExit();
    return 0;
}
