/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    template<Scalar T, int Major0, int Major1, int Major2>
    void gemm(benchmark::State& state) {
        RandomSource::getInstance().reseed();

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
        RandomSource::getInstance().reseed();

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

#ifndef PHYSICA_EmitLLVM
    template<Scalar T, int Major0, int Major1, int Major2>
    void gemm_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

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
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> M0& { expr.assign_base(result); return result; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> M0& { expr.assign_mkl(result); return result; }, "mkl"),
                        backend([&](auto&&...) -> M0& { expr.assign(result); return result; }, "dispatch"));
    }
#endif
}

using enum MatrixMajor::Option;
BENCHMARK(gemm<float64, Col, Col, Col>)->Name("GEMM CCC base")
    ->Args({2, 2, 2})
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});
BENCHMARK(gemm<float64, Col, Row, Col>)->Name("GEMM CRC base")
    ->Args({2, 2, 2})
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});
BENCHMARK(gemm<float64, Col, Col, Row>)->Name("GEMM CCR base")
    ->Args({2, 2, 2})
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});
BENCHMARK(gemm<float64, Col, Row, Row>)->Name("GEMM CRR base")
    ->Args({2, 2, 2})
    ->Args({4, 4, 4})
    ->Args({8, 8, 8})
    ->Args({16, 16, 16})
    ->Args({64, 64, 64})
    ->Args({256, 256, 256})
    ->Args({1024, 1024, 1024});

BENCHMARK(gemm_trans)->Name("GEMM trans");

#ifndef PHYSICA_EmitLLVM
BENCHMARK(gemm_compare<float64, Col, Col, Col>)->Name("GEMM CCC")->UseManualTime()->Args({2, 2, 2})->Args({4, 4, 4})->Args({8, 8, 8})->Args({16, 16, 16})->Args({64, 64, 64})->Args({256, 256, 256})->Args({1024, 1024, 1024});
BENCHMARK(gemm_compare<float64, Col, Row, Col>)->Name("GEMM CRC")->UseManualTime()->Args({2, 2, 2})->Args({4, 4, 4})->Args({8, 8, 8})->Args({16, 16, 16})->Args({64, 64, 64})->Args({256, 256, 256})->Args({1024, 1024, 1024});
BENCHMARK(gemm_compare<float64, Col, Col, Row>)->Name("GEMM CCR")->UseManualTime()->Args({2, 2, 2})->Args({4, 4, 4})->Args({8, 8, 8})->Args({16, 16, 16})->Args({64, 64, 64})->Args({256, 256, 256})->Args({1024, 1024, 1024});
BENCHMARK(gemm_compare<float64, Col, Row, Row>)->Name("GEMM CRR")->UseManualTime()->Args({2, 2, 2})->Args({4, 4, 4})->Args({8, 8, 8})->Args({16, 16, 16})->Args({64, 64, 64})->Args({256, 256, 256})->Args({1024, 1024, 1024});
#endif
