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
    void asum(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const int64_t size = state.range(0);
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            PHYSICA_BENCH(auto y = x.norm1_base());
            benchmark::DoNotOptimize(x);
            benchmark::DoNotOptimize(y);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void asum_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const int64_t size = state.range(0);
        auto x = VectorND<T>::template random_uniform<RandomSource>(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) { return x.norm1_base(); }, "base"),
                        backend<HasMKL()>([&](auto&&...) { return x.norm1_mkl(); }, "mkl"),
                        backend([&](auto&&...) { return x.norm1(); }, "dispatch"));
    }
#endif
}

BENCHMARK(asum)->Name("asum base")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(16384);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(asum_compare<float64>)->Name("asum")->UseManualTime()->Arg(2)->Arg(4)->Arg(8)->Arg(16)->Arg(64)->Arg(256)->Arg(1024)->Arg(16384);
#endif
