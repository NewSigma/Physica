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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    void assign(benchmark::State& state) {
        using T = float64;
        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        MatrixND<T> buffer(order);
        for (auto _ : state) {
            m.tril_unit().assign(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void assign_base(benchmark::State& state) {
        using T = float64;
        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        MatrixND<T> buffer(order);
        for (auto _ : state) {
            m.tril_unit().assign_base(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(assign)->Name("Trig assign")
    ->Arg(256)
    ->Arg(1024)
    ->Arg(8192);

BENCHMARK(assign_base)->Name("Trig assign base")
    ->Arg(256)
    ->Arg(1024)
    ->Arg(8192);
