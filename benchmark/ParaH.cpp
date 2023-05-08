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
#include "Physica/Core/Physics/MD/ForceModel/PairModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Physics/MD/KineticModel/PeriodicModel.h"
#include "Physica/Utils/Random.h"
#include "Physica/Utils/BenchmarkHelper.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;
using PosScalarType = ScalarType;
using ThermostatType = DoubleThermo<ScalarType, PosScalarType>;
using KineticModel = PeriodicModel<ScalarType, PosScalarType, 3>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

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

ScalarType force(ScalarType r) {
    constexpr double alpha = 1.713;
    constexpr double beta = 1.5671;
    constexpr double gamma = 0.00993;
    constexpr double cutoff = 8.32;
    constexpr double c6 = 12.14;
    constexpr double c8 = 215.2;
    constexpr double c9 = 143.1;
    constexpr double c10 = 4813.9;

    const ScalarType r2 = square(r);
    ScalarType result = exp(-r2 * gamma - r * beta + alpha) * (r * (gamma * 2) + beta);
    const ScalarType rep_r = reciprocal(r);
    const ScalarType rep_r2 = square(rep_r);
    const ScalarType rep_r4 = square(rep_r2);
    const ScalarType rep_r6 = rep_r4 * rep_r2;
    const ScalarType rep_r7 = rep_r6 * rep_r;
    const ScalarType rep_r9 = rep_r7 * rep_r2;
    const ScalarType rep_r10 = rep_r6 * rep_r4;
    const ScalarType rep_r11 = rep_r7 * rep_r4;

    const ScalarType d_g = rep_r10 * (9 * c9) - rep_r11 * (10 * c10) - rep_r9 * (8 * c8) - rep_r7 * (6 * c6);
    if (r < cutoff) {
        const ScalarType rep_r8 = square(rep_r4);
        const ScalarType g = rep_r6 * c6 + rep_r8 * c8 - rep_r9 * c9 + rep_r10 * c10;
        const ScalarType f_cutoff = exp(-square(rep_r * cutoff - 1));
        result += (d_g + g * rep_r * (rep_r2 * (2 * cutoff * cutoff) - rep_r * (2 * cutoff))) * f_cutoff;
    }
    else {
        result += d_g;
    }
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
    Cycler::init();
    std::mt19937 gen(3438603950906262893);

    ThreadPool::numThreadRequired = 4;
    ThreadPool& pool = ThreadPool::getInstance();
    {
        auto rpmd = makeSystem(gen);
        const ThermostatType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        rpmd.initMomentum(gen);

        PairModel<ScalarType, PosScalarType, decltype(&force)> pair(ScalarType(pair_cutoff), force, pot_functor);
        rpmd.updateForce<decltype(pair), ThreadExecutor>(pair);

        //ProfilerStart("profiler.dat");
        auto timeuse = Benchmark::run([&]() {
            rpmd.nvt_step_for<ThermostatType, decltype(gen), KineticModel, decltype(pair), ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-13),
                thermo,
                gen,
                kineticModel,
                pair);
        }, 8, 20);
        //ProfilerStop();
        std::cout << "4 Threads time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    pool.shouldExit();
    return 0;
}
