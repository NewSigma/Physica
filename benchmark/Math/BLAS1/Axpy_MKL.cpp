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
    void axpy_mkl(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            (x * T(2)).assign_add_mkl(buffer);
            benchmark::DoNotOptimize(x);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(2);
BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(4);
BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(8);
BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(16);
BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(64);
BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(256);
BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(1024);
BENCHMARK(axpy_mkl)->Name("axpy mkl")->Arg(16384);