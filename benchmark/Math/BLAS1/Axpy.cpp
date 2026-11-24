/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void axpy(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = x * T(2);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_add_base(buffer));
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void axpy_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        auto x = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = x * T(2);
        VectorND<T> buffer(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign_add_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { expr.assign_add_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign_add(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(axpy)->Name("axpy base")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(16384);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(axpy_compare<float64>)->Name("axpy")->UseManualTime()->Arg(2)->Arg(4)->Arg(8)->Arg(16)->Arg(64)->Arg(256)->Arg(1024)->Arg(16384);
#endif
