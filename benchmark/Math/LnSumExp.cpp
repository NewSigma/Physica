/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void lnSumExp(benchmark::State& state) {
        // Benchmark we create a temporary for lnSumExp calculation
        // This pattern is extracted from classical ising model
        RandomSource::getInstance().reseed();

        using T = float32;
        const int64_t size = state.range(0);
        const T factor = 1.2;
        const auto mat = MatrixND<T>::random_uniform<RandomSource>(size, size);
        const auto vecA = VectorND<T>::random_uniform<RandomSource>(size);
        const auto vecB = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = factor * hadamard(vecA, square(mat.transpose() * vecB));
        for (auto _ : state) {
            PHYSICA_BENCH(auto x = expr.lnSumExp());
            benchmark::DoNotOptimize(x);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(lnSumExp)->Name("lnSumExp")->Arg(512);
