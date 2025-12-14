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
using RandomSource = Random<MT19937>;

namespace {
    template<Scalar T>
    void add(benchmark::State& state) {
        const int64_t size = CacheSizes[state.range(0)] / sizeof(T);
        const VectorND<T> a = VectorND<T>::template random_uniform<RandomSource>(size);
        const VectorND<T> b = VectorND<T>::template random_uniform<RandomSource>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            (a + b).assign_mkl(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(add<float32>)->Name("add float32 mkl")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(add<float64>)->Name("add float64 mkl")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
