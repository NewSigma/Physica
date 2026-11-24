/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"
#include "Benchmark.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    template<bool Left>
    void kernel(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t size = state.range(0);
        auto mat = DenseMatrix<T>::template random_uniform<RandomSource>(size, size);
        const auto vec = VectorND<T>::template random_uniform<RandomSource>(size);
        for (auto _ : state) {
            if constexpr (Left)
                PHYSICA_BENCH(applyHouseholder(T(2), vec, mat));
            else
                PHYSICA_BENCH(applyHouseholder(mat, T(2), vec));
            benchmark::DoNotOptimize(mat);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(kernel<true>)->Name("applyHouseHolder left")->Arg(512);
BENCHMARK(kernel<false>)->Name("applyHouseHolder right")->Arg(512);
