/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/FreqDQMC.h"
#include "Benchmark.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;
using RandomSource = Random<MCG, std::mt19937::default_seed>;
constexpr int Dim = 2;
constexpr T HoppingT = 1;
constexpr T RepelU = 4;
constexpr T Beta = 8;
constexpr int NumSiteX = 4;
constexpr int NumSiteY = 4;
constexpr int NumSplit = Beta.toMachine() * 8;

namespace {
    void kernel(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, NumSplit);
        auto dqmc = FreqDQMC<Tc>(params, 2);
        dqmc.step_random<RandomSource>();
        for (auto _ : state) {
            PHYSICA_BENCH(dqmc.step<RandomSource>(););
            benchmark::DoNotOptimize(dqmc);
        }
    }
}

BENCHMARK(kernel)->Name("FreqDQMC");
