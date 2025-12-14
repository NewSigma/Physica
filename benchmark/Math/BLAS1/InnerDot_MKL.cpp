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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

namespace {
    template<Scalar T>
    void innerDot_mkl(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        auto dot = InnerDot(v1, v2);
        for (auto _ : state) {
            auto y = dot.calc_mkl();
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }
}

BENCHMARK(innerDot_mkl<float32>)->Name("innerDot mkl float32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_mkl<float64>)->Name("innerDot mkl float64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_mkl<cfloat32>)->Name("innerDot mkl cfloat32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(innerDot_mkl<cfloat64>)->Name("innerDot mkl cfloat64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);
