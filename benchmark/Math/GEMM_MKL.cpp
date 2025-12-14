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
using RandomSource = Random<MT19937>;

namespace {
    template<int Order = 0>
    void gemm_mkl(benchmark::State& state) {
        static_assert(Order >= 0);
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixOption::Col, Order, Order>;
        const size_t order = state.range(0);
        const auto m1 = MatrixType::template random_uniform<RandomSource>(order, order);
        const auto m2 = MatrixType::template random_uniform<RandomSource>(order, order);
        MatrixType m(order, order);
        for (auto _ : state) {
            (m1 * m2).assign_mkl(m);
            benchmark::DoNotOptimize(m);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(gemm_mkl<2>)->Name("GEMM mkl")->Arg(2);
BENCHMARK(gemm_mkl<4>)->Name("GEMM mkl")->Arg(4);
BENCHMARK(gemm_mkl<8>)->Name("GEMM mkl")->Arg(8);
BENCHMARK(gemm_mkl<16>)->Name("GEMM mkl")->Arg(16);
BENCHMARK(gemm_mkl<32>)->Name("GEMM mkl")->Arg(32);
BENCHMARK(gemm_mkl)->Name("GEMM mkl")
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024);
