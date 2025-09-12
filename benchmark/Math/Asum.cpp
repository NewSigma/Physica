/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

namespace {
    void asum(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(x.norm1());
            benchmark::ClobberMemory();
        }
    }

    void asum_base(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(x.norm1_base());
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(asum)->Name("asum")->Arg(2);
BENCHMARK(asum)->Name("asum")->Arg(4);
BENCHMARK(asum)->Name("asum")->Arg(8);
BENCHMARK(asum)->Name("asum")->Arg(16);
BENCHMARK(asum)->Name("asum")->Arg(64);
BENCHMARK(asum)->Name("asum")->Arg(256);
BENCHMARK(asum)->Name("asum")->Arg(1024);
BENCHMARK(asum)->Name("asum")->Arg(16384);

BENCHMARK(asum_base)->Name("asum base")->Arg(2);
BENCHMARK(asum_base)->Name("asum base")->Arg(4);
BENCHMARK(asum_base)->Name("asum base")->Arg(8);
BENCHMARK(asum_base)->Name("asum base")->Arg(16);
BENCHMARK(asum_base)->Name("asum base")->Arg(64);
BENCHMARK(asum_base)->Name("asum base")->Arg(256);
BENCHMARK(asum_base)->Name("asum base")->Arg(1024);
BENCHMARK(asum_base)->Name("asum base")->Arg(16384);
