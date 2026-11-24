/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void copy(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> y(size);
        for (auto _ : state) {
            PHYSICA_BENCH(x.assign_base(y));
            benchmark::DoNotOptimize(y);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void copy_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        auto x = VectorND<T>::template random_uniform<RandomSource>(size);
        VectorND<T> y(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { x.assign_base(y); return y; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { x.assign_mkl(y); return y; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { x.assign(y); return y; }, "dispatch"));
    }
#endif
}

BENCHMARK(copy)->Name("copy base")
    ->Arg(4)
    ->Arg(8)
    ->Arg(4UL * 1024UL / sizeof(float64))
    ->Arg(64UL * 1024UL / sizeof(float64))
    ->Arg(1024UL * 1024UL / sizeof(float64))
    ->Arg(16UL * 1024UL * 1024UL / sizeof(float64));

#ifndef PHYSICA_EmitLLVM
BENCHMARK(copy_compare<float64>)->Name("copy")->UseManualTime()->Arg(4)->Arg(8)->Arg(4UL * 1024UL / sizeof(float64))->Arg(64UL * 1024UL / sizeof(float64))->Arg(1024UL * 1024UL / sizeof(float64))->Arg(16UL * 1024UL * 1024UL / sizeof(float64));
#endif
