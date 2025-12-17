/*
 * Copyright 2025 Weibo He.
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
using RandomSource = Random<MCG>;

namespace {
    template<Scalar T>
    void square(benchmark::State& state) {
        const int64_t size = CacheSizes[state.range(0)] / sizeof(T);
        const VectorND<T> x = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = square(x);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            [[clang::noinline]] expr.assign(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    template<Scalar T>
    void square_base(benchmark::State& state) {
        const int64_t size = CacheSizes[state.range(0)] / sizeof(T);
        const VectorND<T> x = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = square(x);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            [[clang::noinline]] expr.assign_base(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(square<float32>)->Name("square float32")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(square<float64>)->Name("square float64")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);

BENCHMARK(square_base<float32>)->Name("square float32 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(square_base<float64>)->Name("square float64 base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
