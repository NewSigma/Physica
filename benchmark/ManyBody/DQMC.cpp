/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/DQMC.h"
#include "Benchmark.h"

using namespace Physica;
using T = float64;
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
        auto dqmc = DQMC<T>(params);
        dqmc.step_random<RandomSource>();
        for (auto _ : state) {
            PHYSICA_BENCH(dqmc.step<RandomSource>(););
            benchmark::DoNotOptimize(dqmc);
        }
    }
}

BENCHMARK(kernel)->Name("DQMC");
