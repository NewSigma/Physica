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
#include <benchmark/benchmark.h>
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Physics/MD/ForceModel/CPUGPUModel.cuh"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using Physica::Dynamic;
using ScalarType = float32;
using MDType = RPMD<ScalarType, 3, Physica::Dynamic, PageLockedAllocator<ScalarType>>;
using MDCellType = MDType::MDCellType;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using RandomSource = Random<MT19937, 10000>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

static MDType makeSystem(size_t numMolecular) {
    MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    auto pos = MDCellType::PositionMatrix::random_uniform<RandomSource>(numMolecular, 3);
    MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);
    return MDType(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}

static void bench(benchmark::State& state) {
    ThreadPool::numThreadRequired = 8;
    KineticModel kineticModel(temperatureT, numReplica);
    using HostModel = SilveraGoldman<ScalarType, true, false>;

    using DeviceModel = SilveraGoldman<ScalarType, true, true>;
    using ForceModel = CPUGPUModel<HostModel, DeviceModel>;
    const size_t numMolecular = state.range();
    MDType rpmd = makeSystem(numMolecular);
    rpmd.initMomentum<KineticModel, RandomSource>();
    ForceModel forceModel(4, HostModel(pair_cutoff), numMolecular, pair_cutoff);
    for (auto _ : state)
        rpmd.nve_step<KineticModel, ForceModel, AutoExecutor>(kineticModel, forceModel);
}

BENCHMARK(bench)->Name("ParaH auto")->Unit(benchmark::kMillisecond)->Arg(108)->Arg(256)->Arg(500)->Arg(864);
