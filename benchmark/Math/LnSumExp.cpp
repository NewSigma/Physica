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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    void lnSumExp(benchmark::State& state) {
        // Benchmark we create a temporary for lnSumExp calculation
        // This pattern is extracted from classical ising model
        using T = float32;
        const int64_t size = state.range(0);
        const T factor = 1.2;
        const auto mat = MatrixND<T>::random_uniform<RandomSource>(size, size);
        const auto vecA = VectorND<T>::random_uniform<RandomSource>(size);
        const auto vecB = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = factor * hadamard(vecA, square(mat.transpose() * vecB));
        for (auto _ : state) {
            PHYSICA_BENCH(auto x = expr.lnSumExp());
            benchmark::DoNotOptimize(x);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(lnSumExp)->Name("lnSumExp")->Arg(512);
