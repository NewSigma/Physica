/*
 * Copyright 2024-2026 Weibo He.
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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    template<Scalar T>
    void innerDot(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        auto dot = InnerDot(v1, v2);
        for (auto _ : state) {
            T y = dot.calc();
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }

    template<Scalar T>
    void innerDot_base(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        auto dot = InnerDot(v1, v2);
        for (auto _ : state) {
            T y{};
            PHYSICA_BENCH(y = dot.calc_base());
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }
    // Benchmark that we correct lower complex-real vector inner dot
    void innerDot_complex_real(benchmark::State& state) {
        const VectorND<cfloat64> v1 = VectorND<cfloat64>::template random_uniform<RandomSource>(1024);
        const VectorND<float64> v2 = VectorND<float64>::template random_uniform<RandomSource>(1024);
        auto dot = InnerDot(v1, v2);
        for (auto _ : state) {
            cfloat64 y{};
            PHYSICA_BENCH(y = dot.calc_base());
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }
}

BENCHMARK(innerDot<float32>)->Name("innerDot float32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot<float64>)->Name("innerDot float64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot<cfloat32>)->Name("innerDot cfloat32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot<cfloat64>)->Name("innerDot cfloat64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);
// Baseline
BENCHMARK(innerDot_base<float32>)->Name("innerDot base float32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_base<float64>)->Name("innerDot base float64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_base<cfloat32>)->Name("innerDot base cfloat32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_base<cfloat64>)->Name("innerDot base cfloat64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_complex_real)->Name("innerDot complex-real");
