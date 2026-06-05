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
    void dot(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        auto dot = Dot(v1, v2);
        for (auto _ : state) {
            T y{};
            PHYSICA_BENCH(y = dot.calc_base());
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }
    // Benchmark that we correct lower complex-real vector dot
    void dot_complex_real(benchmark::State& state) {
        const VectorND<cfloat64> v1 = VectorND<cfloat64>::template random_uniform<RandomSource>(1024);
        const VectorND<float64> v2 = VectorND<float64>::template random_uniform<RandomSource>(1024);
        auto dot = Dot(v1, v2);
        for (auto _ : state) {
            cfloat64 y{};
            PHYSICA_BENCH(y = dot.calc_base());
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }
}

BENCHMARK(dot<float32>)->Name("dot base float32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot<float64>)->Name("dot base float64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot<cfloat32>)->Name("dot base cfloat32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot<cfloat64>)->Name("dot base cfloat64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot_complex_real)->Name("dot complex-real");
