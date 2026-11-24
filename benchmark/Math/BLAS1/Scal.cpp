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
    void scal(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        bool flag = true;
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            auto expr = x * T(flag ? std::numbers::pi : (1.0 / std::numbers::pi));
            PHYSICA_BENCH(expr.assign_base(x));
            benchmark::DoNotOptimize(x);
            benchmark::ClobberMemory();
            flag = !flag;
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void scal_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        auto x = VectorND<T>::template random_uniform<RandomSource>(size);
        bool flag = true;
        const auto factor = [&]() {
            // Dynamic factor ensures we do not overflow or creep towards zero
            const T alpha = flag ? T(std::numbers::pi) : T(1.0 / std::numbers::pi);
            flag = !flag;
            return alpha;
        };
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { (x * factor()).assign_base(x); return x; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { (x * factor()).assign_mkl(x); return x; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { (x * factor()).assign(x); return x; }, "dispatch"));
    }
#endif
}

BENCHMARK(scal)->Name("scal base")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(16384);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(scal_compare<float64>)->Name("scal")->UseManualTime()->Arg(2)->Arg(4)->Arg(8)->Arg(16)->Arg(64)->Arg(256)->Arg(1024)->Arg(16384);
#endif
