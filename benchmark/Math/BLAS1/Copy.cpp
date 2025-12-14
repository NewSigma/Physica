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
    void copy(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        const auto x = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> y(size);
        for (auto _ : state) {
            [[clang::noinline]] x.assign(y);
            benchmark::DoNotOptimize(y);
            benchmark::ClobberMemory();
        }
    }

    void copy_base(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> y(size);
        for (auto _ : state) {
            [[clang::noinline]] x.assign_base(y);
            benchmark::DoNotOptimize(y);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(copy)->Name("copy")
    ->Arg(4)
    ->Arg(8)
    ->Arg(4UL * 1024UL / sizeof(float64))
    ->Arg(64UL * 1024UL / sizeof(float64))
    ->Arg(1024UL * 1024UL / sizeof(float64))
    ->Arg(16UL * 1024UL * 1024UL / sizeof(float64));

BENCHMARK(copy_base)->Name("copy base")
    ->Arg(4)
    ->Arg(8)
    ->Arg(4UL * 1024UL / sizeof(float64))
    ->Arg(64UL * 1024UL / sizeof(float64))
    ->Arg(1024UL * 1024UL / sizeof(float64))
    ->Arg(16UL * 1024UL * 1024UL / sizeof(float64));
