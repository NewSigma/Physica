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
#include "Exp.h"

using namespace Physica;

namespace {
    template<Scalar T>
    void exp(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> data = makeData<T>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            exp(data).assign_mkl(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(exp<float32>)->Name("exp float32 mkl")->Unit(benchmark::kNanosecond)->Arg(8);
BENCHMARK(exp<float32>)->Name("exp float32 mkl")->Unit(benchmark::kNanosecond)->Arg(64);
BENCHMARK(exp<float32>)->Name("exp float32 mkl")->Unit(benchmark::kNanosecond)->Arg(512);
BENCHMARK(exp<float32>)->Name("exp float32 mkl")->Unit(benchmark::kNanosecond)->Arg(1024);
BENCHMARK(exp<float64>)->Name("exp float64 mkl")->Unit(benchmark::kNanosecond)->Arg(8);
BENCHMARK(exp<float64>)->Name("exp float64 mkl")->Unit(benchmark::kNanosecond)->Arg(16);
BENCHMARK(exp<float64>)->Name("exp float64 mkl")->Unit(benchmark::kNanosecond)->Arg(32);
BENCHMARK(exp<float64>)->Name("exp float64 mkl")->Unit(benchmark::kNanosecond)->Arg(64);
