/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

namespace {
    template<Scalar T>
    void innerDot(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        const auto dot = InnerDot(v1, v2);
        for (auto _ : state)
            benchmark::DoNotOptimize(dot.calc());
    }

    template<Scalar T>
    void innerDot_base(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        const auto dot = InnerDot(v1, v2);
        for (auto _ : state)
            benchmark::DoNotOptimize(dot.calc_base());
    }
}

BENCHMARK(innerDot<float32>)->Name("innerDot float32")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot<float64>)->Name("innerDot float64")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot<cfloat32>)->Name("innerDot cfloat32")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot<cfloat64>)->Name("innerDot cfloat64")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);
// Baseline
BENCHMARK(innerDot_base<float32>)->Name("innerDot base float32")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_base<float64>)->Name("innerDot base float64")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_base<cfloat32>)->Name("innerDot base cfloat32")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_base<cfloat64>)->Name("innerDot base cfloat64")->Unit(benchmark::kNanosecond)
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);
