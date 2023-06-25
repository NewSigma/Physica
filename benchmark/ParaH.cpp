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
#include <gperftools/profiler.h>
#include <iostream>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Utils/Random.h"
#include "Physica/Utils/BenchmarkHelper.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;
using PosScalarType = ScalarType;
using ThermostatType = DoubleThermo<ScalarType, PosScalarType>;
using KineticModel = FreeModel<ScalarType, PosScalarType, 3>;
using ForceModel = SilveraGoldman<ScalarType, PosScalarType>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

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
    std::mt19937 gen(3438603950906262893);

    ThreadPool::numThreadRequired = 4;
    ThreadPool& pool = ThreadPool::getInstance();
    {
        auto rpmd = makeSystem(gen);
        const ThermostatType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        rpmd.initMomentum(gen);

        ForceModel forceModel(pair_cutoff);
        rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);

        //ProfilerStart("profiler.dat");
        auto timeuse = Benchmark::run([&]() {
            rpmd.nvt_step_for<ThermostatType, decltype(gen), KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-13),
                thermo,
                gen,
                kineticModel,
                forceModel);
        }, 8, 20);
        //ProfilerStop();
        std::cout << "4 Threads time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    pool.shouldExit();
    return 0;
}
