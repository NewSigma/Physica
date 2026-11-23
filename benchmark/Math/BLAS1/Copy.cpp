/*
 * Copyright 2025-2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void copy(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> y(size);
        for (auto _ : state) {
            PHYSICA_BENCH(x.assign_base(y));
            benchmark::DoNotOptimize(y);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void copy_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        auto x = VectorND<T>::template random_uniform<RandomSource>(size);
        VectorND<T> y(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { x.assign_base(y); return y; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { x.assign_mkl(y); return y; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { x.assign(y); return y; }, "dispatch"));
    }
#endif
}

BENCHMARK(copy)->Name("copy base")
    ->Arg(4)
    ->Arg(8)
    ->Arg(4UL * 1024UL / sizeof(float64))
    ->Arg(64UL * 1024UL / sizeof(float64))
    ->Arg(1024UL * 1024UL / sizeof(float64))
    ->Arg(16UL * 1024UL * 1024UL / sizeof(float64));

#ifndef PHYSICA_EmitLLVM
BENCHMARK(copy_compare<float64>)->Name("copy")->UseManualTime()->Arg(4)->Arg(8)->Arg(4UL * 1024UL / sizeof(float64))->Arg(64UL * 1024UL / sizeof(float64))->Arg(1024UL * 1024UL / sizeof(float64))->Arg(16UL * 1024UL * 1024UL / sizeof(float64));
#endif
