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
    template<Scalar T>
    void abs(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeVectorSize(state.range(0), sizeof(T));
        const VectorND<T> x = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = abs(x);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_base(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void abs_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeVectorSize(state.range(0), sizeof(T));
        const VectorND<T> x = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = abs(x);
        VectorND<T> buffer(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { expr.assign_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(abs<float32>)->Name("abs float32 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(abs<float64>)->Name("abs float64 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(abs_compare<float32>)->Name("abs float32")->UseManualTime()->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(abs_compare<float64>)->Name("abs float64")->UseManualTime()->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
#endif
