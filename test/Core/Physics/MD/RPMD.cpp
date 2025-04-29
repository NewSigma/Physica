/*
 * Copyright 2022-2024 Weibo He.
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

using namespace Physica;
using RandomSource = Random<MT19937, 3438603950906262893>;
using ScalarType = float64;
using KineticModel = FreeModel<ScalarType, 3, Physica::Dynamic, RPMDIntegrator::Exact>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

RPMD<ScalarType> makeSystem() {
    using MDCellType = RPMD<ScalarType>::MDCellType;
    MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    auto pos = MDCellType::PositionMatrix::random_uniform<RandomSource>(numMolecular, 3);
    MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}

bool testDriftMomentum(double precision) {
    auto rpmd = makeSystem();
    rpmd.initMomentum<KineticModel, RandomSource>();
    for (int i = 0; i < 3; ++i) {
        ScalarType sum = 0;
        for (size_t j = i; j < rpmd.getDOF(); j += 3)
            sum += rpmd.getPhaseMatrix().row(j).sum();
        if (!scalarNear(sum, ScalarType(0), precision))
            return false;
    }
    return true;
}

bool testCalcKinetic(double precision) {
    using ForceModel = SilveraGoldman<ScalarType, true, false>;

    auto rpmd = makeSystem();
    rpmd.initMomentum<KineticModel, RandomSource>();
    ForceModel forceModel(pair_cutoff);
    rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);

    const ScalarType kinetic1 = rpmd.calcKinetic<KineticModel>();
    ScalarType kinetic2 = 0;
    for (size_t i = 0; i < rpmd.getDOF(); ++i)
        kinetic2 += rpmd.calcKinetic<KineticModel>(i);
    return scalarNear(kinetic1, kinetic2, precision);
}
/**
 * Reference:
 * [1] J. Chem. Phys. 122, 184503 (2005); https://doi.org/10.1063/1.1893956
 */
void testMDRun() {
    using ForceModel = SilveraGoldman<ScalarType, true, false>;
    using ThermoType = DoubleThermo<KineticModel>;
    ScalarType mean = 0;
    ScalarType var = 0;
    {
        const ThermoType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        ForceModel forceModel(pair_cutoff);

        auto rpmd = makeSystem();
        rpmd.initMomentum<KineticModel, RandomSource>();
        for (unsigned int i = 0; i < 6; ++i) {
            ScalarType temp = 0;
            rpmd.nvt_step_for<ThermoType, RandomSource, KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-12), thermo, kineticModel, forceModel);

            for (unsigned int j = 0; j < 100; ++j) {
                rpmd.nvt_step<ThermoType, RandomSource, KineticModel, ForceModel, ThreadExecutor>(thermo, kineticModel, forceModel);
                toNextMean(temp, j, rpmd.calcKinetic<KineticModel>());
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
