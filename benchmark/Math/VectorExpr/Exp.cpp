/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Exp.h"
#include "Benchmark.h"

using namespace Physica;

namespace {
    template<Scalar T>
    void exp(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> data = makeData<T>(size);
        auto expr = exp(data);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_base(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void exp_compare(benchmark::State& state) {
        const int64_t size = state.range(0);
        const VectorND<T> data = makeData<T>(size);
        auto expr = exp(data);
        VectorND<T> buffer(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { expr.assign_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(exp<float32>)->Name("exp float32 base")
    ->Arg(8)
    ->Arg(64)
    ->Arg(512)
    ->Arg(1024);
BENCHMARK(exp<float64>)->Name("exp float64 base")
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(exp_compare<float32>)->Name("exp float32")->UseManualTime()->Arg(8)->Arg(64)->Arg(512)->Arg(1024);
BENCHMARK(exp_compare<float64>)->Name("exp float64")->UseManualTime()->Arg(8)->Arg(16)->Arg(32)->Arg(64);
#endif
