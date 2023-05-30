/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Physics/MD/KineticModel/PeriodicModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using PosScalarType = ScalarType;
using ThermostatType = DoubleThermo<ScalarType, PosScalarType>;
using KineticModel = PeriodicModel<ScalarType, PosScalarType, 3>;
using ForceModel = SilveraGoldman<ScalarType, PosScalarType>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

bool testDriftMomentum(const RPMD<ScalarType, PosScalarType>& rpmd, double precision) {
    for (int i = 0; i < 3; ++i) {
        ScalarType sum = 0;
        for (size_t j = i; j < rpmd.getDOF(); j += 3)
            sum += rpmd.getPhaseMatrix().row(j).asVector().sum();
        if (!scalarNear(sum, ScalarType::Zero(), precision))
            return false;
    }
    return true;
}

ScalarType pot_functor(ScalarType r) {
    constexpr double alpha = 1.713;
    constexpr double beta = 1.5671;
    constexpr double gamma = 0.00993;
    constexpr double cutoff = 8.32;
    constexpr double c6 = 12.14;
    constexpr double c8 = 215.2;
    constexpr double c9 = 143.1;
    constexpr double c10 = 4813.9;

    const ScalarType r2 = square(r);
    ScalarType result = exp(-r2 * gamma - r * beta + alpha);
    const ScalarType rep_r = reciprocal(r);
    const ScalarType rep_r2 = square(rep_r);
    const ScalarType rep_r4 = square(rep_r2);
    const ScalarType rep_r6 = rep_r4 * rep_r2;
    const ScalarType rep_r8 = square(rep_r4);
    const ScalarType rep_r9 = rep_r8 * rep_r;
    const ScalarType rep_r10 = rep_r6 * rep_r4;
    const ScalarType g = rep_r6 * c6 + rep_r8 * c8 - rep_r9 * c9 + rep_r10 * c10;

    if (r < cutoff) {
        const ScalarType f_cutoff = exp(-square(rep_r * cutoff - 1));
        result -= g * f_cutoff;
    }
    else
        result -= g;
    return result;
}

template<class RandomGenerator>
RPMD<ScalarType, PosScalarType> makeSystem(RandomGenerator& gen) {
    using MDCellType = typename RPMD<ScalarType, PosScalarType>::MDCellType;
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType, PosScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
int main() {
    ScalarType mean = 0;
    ScalarType var = 0;

    ThreadPool::numThreadRequired = 4;
    ThreadPool& pool = ThreadPool::getInstance();
    {
        std::mt19937 gen(3438603950906262893);
        auto rpmd = makeSystem(gen);
        const ThermostatType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        rpmd.initMomentum(gen);
        if (!testDriftMomentum(rpmd, 1E-12))
            return 1;

        ForceModel forceModel(pair_cutoff);
        rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);

        for (unsigned int i = 0; i < 6; ++i) {
            ScalarType temp = 0;
            rpmd.nvt_step_for<ThermostatType, decltype(gen), KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-12),
                thermo,
                gen,
                kineticModel,
                forceModel);

            for (unsigned int j = 0; j < 100; ++j) {
                rpmd.nvt_step<ThermostatType, decltype(gen), KineticModel, ForceModel, ThreadExecutor>(thermo, gen, kineticModel, forceModel);
                toNextMean(temp, j, rpmd.calcKinetic());
            }
            toNextVariance(var, mean, i, temp);
        }
    }
    pool.shouldExit();

    constexpr double answer = 61.8;
    const ScalarType energyPerMol = PhyConst<AU>::temperatureToK(double(mean) / numMolecular);
    if (!scalarNear(energyPerMol, ScalarType(answer), 0.1))
        return 1;
    const ScalarType delta = abs(energyPerMol - answer);
    const ScalarType deviation = PhyConst<AU>::temperatureToK(std::sqrt(double(var))) / numMolecular;
    if (delta > deviation)
        return 1;
    return 0;
}
