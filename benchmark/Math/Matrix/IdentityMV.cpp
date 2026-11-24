/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void kernel1(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        const auto v1 = VectorND<T>::random_uniform<RandomSource>(size);
        const auto v2 = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = v1 - v2 * T(2.1);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            expr.assign(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void kernel2(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        const auto v1 = VectorND<T>::random_uniform<RandomSource>(size);
        const auto v2 = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = v1 - (IdentityMatrix<T>(size) * T(2.1)) * v2;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void kernel3(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        const auto v = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = (IdentityMatrix<T>(size) * T(2.1)) * v;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            expr.assign_add(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(kernel1)->Name("IdentityMV 1")->Arg(256);
BENCHMARK(kernel2)->Name("IdentityMV 2")->Arg(256);
BENCHMARK(kernel3)->Name("IdentityMV 3")->Arg(256);
