/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

namespace {
    void kernel(benchmark::State& state) {
        using T = float64;
        using MatrixType = DenseSymmMatrix<T>;
        const int size = state.range(0);
        const auto m = MatrixType::random_uniform<RandomSource>(size);
        const auto v = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            buffer = m * v;
            benchmark::DoNotOptimize(buffer);
        }
    }
}

BENCHMARK(kernel)->Name("SyMV")->Unit(benchmark::kNanosecond)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(64);
BENCHMARK(kernel)->Name("SyMV")->Unit(benchmark::kMicrosecond)->Arg(256);
BENCHMARK(kernel)->Name("SyMV")->Unit(benchmark::kMicrosecond)->Arg(1024);
BENCHMARK(kernel)->Name("SyMV")->Unit(benchmark::kMillisecond)->Arg(16384);
