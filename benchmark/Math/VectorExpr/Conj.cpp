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
    void conj(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeVectorSize(state.range(0), sizeof(T));
        const VectorND<T> x = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = x.conjugate();
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_base(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void conj_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeVectorSize(state.range(0), sizeof(T));
        const VectorND<T> x = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = x.conjugate();
        VectorND<T> buffer(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { expr.assign_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(conj<cfloat32>)->Name("conj cfloat32 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(conj<cfloat64>)->Name("conj cfloat64 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(conj_compare<cfloat32>)->Name("conj cfloat32")->UseManualTime()->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(conj_compare<cfloat64>)->Name("conj cfloat64")->UseManualTime()->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
#endif
