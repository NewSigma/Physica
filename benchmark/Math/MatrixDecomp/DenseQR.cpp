/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseQR.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void qr(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        DenseQR<T> qr(order, order);
        for (auto _ : state) {
            PHYSICA_BENCH(qr.compute_base(m));
            benchmark::DoNotOptimize(qr);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void qr_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        DenseQR<T> qr(order, order);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> DenseQR<T>& { qr.compute_base(m); return qr; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> DenseQR<T>& { qr.compute_mkl(m); return qr; }, "mkl"),
                        backend([&](auto&&...) -> DenseQR<T>& { qr.compute(m); return qr; }, "dispatch"));
    }
#endif
}

BENCHMARK(qr)->Name("DenseQR base")
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(qr_compare<float64>)->Name("DenseQR")->UseManualTime()->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64)->Arg(256)->Arg(1024);
#endif
