/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Benchmark.h"

using namespace Physica;
using T = float32;
using KineticModel = FreeModel<T, 3, Dynamic, RPMDIntegrator::Exact>;
using RandomSource = Random<MCG, std::mt19937::default_seed>;
using MDType = RPMD<T, 3, Physica::Dynamic, PageLockedAllocator<T>>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

namespace {
    MDType makeSystem(size_t numMolecular) {
        using MDCellType = MDType::MDCellType;
        MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::identity(3);
        auto pos = MDCellType::PositionMatrix::random_uniform<RandomSource>(numMolecular, 3);
        MDCellType::MassVector massVec(numMolecular, mass);
        MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

        const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
        cell.scale(factor);
        return MDType(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
    }

    void bench(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using ForceModel = device_obj<SilveraGoldman<T, true, true>>;
        const int numMolecular = state.range(0);
        KineticModel kineticModel(temperatureT, numReplica);
        MDType rpmd = makeSystem(numMolecular);
        rpmd.initMomentum<KineticModel, RandomSource>();
        ForceModel forceModel(numMolecular, pair_cutoff);
        for (auto _ : state)
            PHYSICA_BENCH(rpmd.nve_step<GPU>(kineticModel, forceModel));
    }
}

BENCHMARK(bench)->Name("ParaH cuda")->Arg(108)->Arg(256)->Arg(500)->Arg(864);
