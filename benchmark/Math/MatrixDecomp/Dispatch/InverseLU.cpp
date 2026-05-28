/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    void lu(benchmark::State& state) {
        using T = float64;
        const size_t order = state.range(0);
        auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        auto lu = DenseLU<T, false>(m);
        auto inv = lu.inv();
        for (auto _ : state) {
            inv.assign(m);
            benchmark::DoNotOptimize(m);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(lu)->Name("DenseLU inv dispatch")
    ->Arg(8)
    ->Arg(64)
    ->Arg(1024);
