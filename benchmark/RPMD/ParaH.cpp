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
#include <benchmark/benchmark.h>
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"

using namespace Physica;
using T = float32;
using KineticModel = FreeModel<T, 3, Physica::Dynamic, RPMDIntegrator::Exact>;
using ForceModel = SilveraGoldman<T, true>;
using RandomSource = Random<MCG>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

namespace {
    RPMD<T> makeSystem(size_t numMolecular) {
        using MDCellType = RPMD<T>::MDCellType;
        MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::identity(3);
        auto pos = MDCellType::PositionMatrix::random_uniform<RandomSource>(numMolecular, 3);
        MDCellType::MassVector massVec(numMolecular, mass);
        MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

        const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
        cell.scale(factor);
        return RPMD<T>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
    }

    void bench(benchmark::State& state) {
        ThreadPool::numThreadRequired = 4;
        KineticModel kineticModel(temperatureT, numReplica);
        ForceModel forceModel(pair_cutoff);

        const int numMolecular = state.range(0);
        auto rpmd = makeSystem(numMolecular);
        rpmd.initMomentum<KineticModel, RandomSource>();
        for (auto _ : state)
            [[clang::noinline]] rpmd.nve_step<Thread>(kineticModel, forceModel);
    }
}

BENCHMARK(bench)->Name("ParaH")->Arg(108)->Arg(256)->Arg(500)->Arg(864);
