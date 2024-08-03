/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Parallel/Executor/CudaExecutor.cuh"

using namespace Physica::Core;
using Physica::Dynamic;
using ScalarType = Scalar<Float>;
using ForceModel = device_obj<SilveraGoldman<ScalarType, true>>;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using ThermoType = DoubleThermo<KineticModel>;
using MDType = RPMD<ScalarType, 3, Physica::Dynamic, Physica::Utils::PageLockedAllocator<ScalarType>>;
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

MDType makeSystem(RandomGenerator& gen) {
    using MDCellType = typename MDType::MDCellType;
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return MDType(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
void testMDRun() {
    ScalarType mean = 0;
    ScalarType var = 0;
    {
        const ThermoType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        ForceModel forceModel(numMolecular, pair_cutoff);
        auto& pool = RandomPoolType::getInstance();
        auto& gen = pool.getGen();
        auto rpmd = makeSystem(gen);
        rpmd.initMomentum<KineticModel, decltype(gen)>(gen);

        for (unsigned int i = 0; i < 6; ++i) {
            ScalarType temp = 0;
            rpmd.nvt_step_for<ThermoType, RandomPoolType, KineticModel, ForceModel, CudaExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-12),
                thermo,
                pool,
                kineticModel,
                forceModel);

            for (unsigned int j = 0; j < 100; ++j) {
                rpmd.nvt_step<ThermoType, RandomPoolType, KineticModel, ForceModel, CudaExecutor>(thermo, pool, kineticModel, forceModel);
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
    testMDRun();
    return 0;
}
