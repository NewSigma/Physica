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
    void nrm2_mkl(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(x.norm2_mkl());
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(2);
BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(4);
BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(8);
BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(16);
BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(64);
BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(256);
BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(1024);
BENCHMARK(nrm2_mkl)->Name("nrm2 mkl")->Arg(16384);
