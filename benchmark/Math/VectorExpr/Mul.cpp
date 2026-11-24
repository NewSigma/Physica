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
    template<Scalar T>
    void mul(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeVectorSize(state.range(0), sizeof(T));
        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> b = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = hadamard(a, b);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_base(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void mul_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeVectorSize(state.range(0), sizeof(T));
        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> b = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = hadamard(a, b);
        VectorND<T> buffer(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { expr.assign_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(mul<float32>)->Name("mul float32 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(mul<float64>)->Name("mul float64 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(mul_compare<float32>)->Name("mul float32")->UseManualTime()->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(mul_compare<float64>)->Name("mul float64")->UseManualTime()->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
#endif
