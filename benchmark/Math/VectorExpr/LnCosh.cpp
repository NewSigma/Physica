/*
 * Copyright 2024-2025 Weibo He.
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
#include <benchmark/benchmark.h>
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    template<Scalar T>
    void lncosh(benchmark::State& state) {
        const VectorND<T> data = VectorND<T>::template random_uniform<RandomSource>(1024);
        VectorND<T> buffer(1024);
        for (auto _ : state) {
            buffer = lncosh(data);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(lncosh<float32>)->Name("lncosh float32");
BENCHMARK(lncosh<float64>)->Name("lncosh float64");
BENCHMARK(lncosh<cfloat32>)->Name("lncosh cfloat32");
BENCHMARK(lncosh<cfloat64>)->Name("lncosh cfloat64");
