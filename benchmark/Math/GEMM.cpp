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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica;
using RandomType = Random<MT19937>;

template<int Order = 0>
static void gemm(benchmark::State& state) {
    static_assert(Order >= 0);
    using T = float64;
    using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Order, Order>;
    const size_t order = state.range(0);
    const auto m1 = MatrixType::template random_uniform<RandomType>(order, order);
    const auto m2 = MatrixType::template random_uniform<RandomType>(order, order);
    MatrixType m(order, order);
    for (auto _ : state) {
        m = m1 * m2;
        benchmark::DoNotOptimize(m);
        benchmark::ClobberMemory();
    }
}

template<int Order = 0>
static void gemm_base(benchmark::State& state) {
    static_assert(Order >= 0);
    using T = float64;
    using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Order, Order>;
    const size_t order = state.range(0);
    const auto m1 = MatrixType::template random_uniform<RandomType>(order, order);
    const auto m2 = MatrixType::template random_uniform<RandomType>(order, order);
    MatrixType m(order, order);
    for (auto _ : state) {
        (m1 * m2).assign_base(m);
        benchmark::DoNotOptimize(m);
        benchmark::ClobberMemory();
    }
}

BENCHMARK(gemm<2>)->Name("GEMM")->Unit(benchmark::kNanosecond)->Arg(2);
BENCHMARK(gemm<4>)->Name("GEMM")->Unit(benchmark::kNanosecond)->Arg(4);
BENCHMARK(gemm<8>)->Name("GEMM")->Unit(benchmark::kNanosecond)->Arg(8);
BENCHMARK(gemm<16>)->Name("GEMM")->Unit(benchmark::kNanosecond)->Arg(16);
BENCHMARK(gemm<32>)->Name("GEMM")->Unit(benchmark::kNanosecond)->Arg(32);
BENCHMARK(gemm)->Name("GEMM")->Unit(benchmark::kMicrosecond)->Arg(64);
BENCHMARK(gemm)->Name("GEMM")->Unit(benchmark::kMillisecond)->Arg(256);
BENCHMARK(gemm)->Name("GEMM")->Unit(benchmark::kMillisecond)->Arg(1024);
BENCHMARK(gemm)->Name("GEMM")->Unit(benchmark::kMillisecond)->Arg(4096);

BENCHMARK(gemm_base<2>)->Name("GEMM base")->Unit(benchmark::kNanosecond)->Arg(2);
BENCHMARK(gemm_base<4>)->Name("GEMM base")->Unit(benchmark::kNanosecond)->Arg(4);
BENCHMARK(gemm_base<8>)->Name("GEMM base")->Unit(benchmark::kNanosecond)->Arg(8);
BENCHMARK(gemm_base<16>)->Name("GEMM base")->Unit(benchmark::kNanosecond)->Arg(16);
BENCHMARK(gemm_base<32>)->Name("GEMM base")->Unit(benchmark::kNanosecond)->Arg(32);
BENCHMARK(gemm_base)->Name("GEMM base")->Unit(benchmark::kMicrosecond)->Arg(64);
BENCHMARK(gemm_base)->Name("GEMM base")->Unit(benchmark::kMillisecond)->Arg(256);
BENCHMARK(gemm_base)->Name("GEMM base")->Unit(benchmark::kMillisecond)->Arg(1024);
//BENCHMARK(gemm_base)->Name("GEMM base")->Unit(benchmark::kMillisecond)->Arg(4096); // Out of time limit
