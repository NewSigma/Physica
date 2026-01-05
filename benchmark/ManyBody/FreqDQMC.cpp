/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/FreqDQMC.h"
#include "Benchmark.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;
using RandomSource = Random<MCG>;
constexpr int Dim = 2;
constexpr T HoppingT = 1;
constexpr T RepelU = 4;
constexpr T Beta = 8;
constexpr int NumSiteX = 4;
constexpr int NumSiteY = 4;
constexpr int NumSplit = Beta.toMachine() * 8;

namespace {
    void kernel(benchmark::State& state) {
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
