/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    template<Scalar T>
    void dot(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        auto dot = Dot(v1, v2);
        for (auto _ : state) {
            T y{};
            PHYSICA_BENCH(y = dot.calc_base());
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }
    // Benchmark that we correct lower complex-real vector dot
    void dot_complex_real(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<cfloat64> v1 = VectorND<cfloat64>::template random_uniform<RandomSource>(1024);
        const VectorND<float64> v2 = VectorND<float64>::template random_uniform<RandomSource>(1024);
        auto dot = Dot(v1, v2);
        for (auto _ : state) {
            cfloat64 y{};
            PHYSICA_BENCH(y = dot.calc_base());
            benchmark::DoNotOptimize(y);
            benchmark::DoNotOptimize(dot);
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void dot_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        const VectorND<T> v1 = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> v2 = VectorND<T>::template random_uniform<RandomSource>(size);
        auto op = Dot(v1, v2);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) { return op.calc_base(); }, "base"),
                        backend<HasMKL()>([&](auto&&...) { return op.calc_mkl(); }, "mkl"),
                        backend([&](auto&&...) { return op.calc(); }, "dispatch"));
    }
#endif
}

BENCHMARK(dot<float32>)->Name("dot base float32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot<float64>)->Name("dot base float64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot<cfloat32>)->Name("dot base cfloat32")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot<cfloat64>)->Name("dot base cfloat64")
    ->Arg(16)
    ->Arg(128)
    ->Arg(1024);

BENCHMARK(dot_complex_real)->Name("dot complex-real");

#ifndef PHYSICA_EmitLLVM
BENCHMARK(dot_compare<float32>)->Name("dot float32")->UseManualTime()->Arg(16)->Arg(128)->Arg(1024);

BENCHMARK(dot_compare<float64>)->Name("dot float64")->UseManualTime()->Arg(16)->Arg(128)->Arg(1024);

BENCHMARK(dot_compare<cfloat32>)->Name("dot cfloat32")->UseManualTime()->Arg(16)->Arg(128)->Arg(1024);

BENCHMARK(dot_compare<cfloat64>)->Name("dot cfloat64")->UseManualTime()->Arg(16)->Arg(128)->Arg(1024);
#endif
