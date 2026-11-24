/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Parallel/Algorithm/Thread.h"
#include "Benchmark.h"

using namespace Physica;
using T = float32;
using KineticModel = FreeModel<T, 3, Physica::Dynamic, RPMDIntegrator::Exact>;
using ForceModel = SilveraGoldman<T, true>;
using RandomSource = Random<MCG, std::mt19937::default_seed>;
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

        const double factor = (std::cbrt(double(numMolecular) * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
        cell.scale(factor);
        return RPMD<T>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
    }

    void bench(benchmark::State& state) {
        RandomSource::getInstance().reseed();
        KineticModel kineticModel(temperatureT, numReplica);
        ForceModel forceModel(pair_cutoff);

        auto rpmd = makeSystem(state.range(0));
        rpmd.initMomentum<KineticModel, RandomSource>();
        for (auto _ : state)
            PHYSICA_BENCH(rpmd.nve_step<Thread>(kineticModel, forceModel));
    }
}

BENCHMARK(bench)->Name("ParaH")->Arg(108)->Arg(256)->Arg(500)->Arg(864);
