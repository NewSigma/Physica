/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        DenseLU<T, false> lu(order);
        for (auto _ : state) {
            PHYSICA_BENCH(lu.compute_base(m));
            benchmark::DoNotOptimize(lu);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void lu_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        DenseLU<T, false> lu(order);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> DenseLU<T, false>& { lu.compute_base(m); return lu; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> DenseLU<T, false>& { lu.compute_mkl(m); return lu; }, "mkl"),
                        backend([&](auto&&...) -> DenseLU<T, false>& { lu.compute(m); return lu; }, "dispatch"));
    }
#endif
}

BENCHMARK(lu)->Name("DenseLU base")
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(lu_compare<float64>)->Name("DenseLU")->UseManualTime()->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64)->Arg(256)->Arg(1024);
#endif
