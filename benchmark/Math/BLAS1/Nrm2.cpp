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
    void nrm2(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            T y{};
            PHYSICA_BENCH(y = x.norm2_base());
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(y);
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void nrm2_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        auto x = VectorND<T>::template random_uniform<RandomSource>(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) { return x.norm2_base(); }, "base"),
                        backend<HasMKL()>([&](auto&&...) { return x.norm2_mkl(); }, "mkl"),
                        backend([&](auto&&...) { return x.norm2(); }, "dispatch"));
    }
#endif
}

BENCHMARK(nrm2)->Name("nrm2 base")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(16384);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(nrm2_compare<float64>)->Name("nrm2")->UseManualTime()->Arg(2)->Arg(4)->Arg(8)->Arg(16)->Arg(64)->Arg(256)->Arg(1024)->Arg(16384);
#endif
