/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void lu(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const size_t order = state.range(0);
        auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        auto lu = DenseLU<T, false>(m);
        auto inv = lu.inv();
        for (auto _ : state) {
            PHYSICA_BENCH(inv.assign_base(m));
            benchmark::DoNotOptimize(m);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void inv_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t order = state.range(0);
        auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        auto lu = DenseLU<T, false>(m);
        auto inv = lu.inv();
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> MatrixND<T>& { inv.assign_base(m); return m; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> MatrixND<T>& { inv.assign_mkl(m); return m; }, "mkl"),
                        backend([&](auto&&...) -> MatrixND<T>& { inv.assign(m); return m; }, "dispatch"));
    }
#endif
}

BENCHMARK(lu)->Name("DenseLU inv base")
    ->Arg(8)
    ->Arg(64)
    ->Arg(1024);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(inv_compare<float64>)->Name("DenseLU inv")->UseManualTime()->Arg(8)->Arg(64)->Arg(1024);
#endif
