/*
 * Copyright 2025-2026 Weibo He.
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
    template<Scalar T, int Major0, int Major1, int Major2>
    void gemm(benchmark::State& state) {
        using M0 = DenseMatrix<T, Major0>;
        using M1 = DenseMatrix<T, Major1>;
        using M2 = DenseMatrix<T, Major2>;

        const size_t m = state.range(0);
        const size_t k = state.range(0);
        const size_t n = state.range(0);
        const auto m1 = M1::template random_uniform<RandomSource>(m, k);
        const auto m2 = M2::template random_uniform<RandomSource>(k, n);
        auto expr = m1 * m2;
        M0 result(m, n);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_base(result));
            benchmark::DoNotOptimize(result);
            benchmark::ClobberMemory();
        }
    }

    void gemm_trans(benchmark::State& state) {
        constexpr int N = 1024;
        const auto m1 = MatrixND<float64>::template random_uniform<RandomSource>(N, N);
        const auto m2 = MatrixND<float64>::template random_uniform<RandomSource>(N, N);
        auto expr = m1 * m2.transpose();
        MatrixND<float64> result(N, N);
        for (auto _ : state) {
            expr.assign(result);
            benchmark::DoNotOptimize(result);
            benchmark::ClobberMemory();
        }
    }
}

using enum MatrixMajor::Option;
BENCHMARK(gemm<float64, Col, Col, Col>)->Name("GEMM CCC base")
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});
BENCHMARK(gemm<float64, Col, Row, Col>)->Name("GEMM CRC base")
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});
BENCHMARK(gemm<float64, Col, Col, Row>)->Name("GEMM CCR base")
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});
BENCHMARK(gemm<float64, Col, Row, Row>)->Name("GEMM CRR base")
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});

BENCHMARK(gemm_trans)->Name("GEMM trans");
