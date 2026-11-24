/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void kernel(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        using MatrixType = DenseSymmMatrix<T>;
        const auto size = state.range(0);
        const auto m = MatrixType::random_uniform<RandomSource>(size);
        const auto v = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = m * v;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
        }
    }
}

BENCHMARK(kernel)->Name("SyMV")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(16384);
