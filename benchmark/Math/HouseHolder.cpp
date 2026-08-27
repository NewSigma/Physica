/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"
#include "Benchmark.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<MCG>;

namespace {
    template<bool Left>
    void kernel(benchmark::State& state) {
        const size_t size = state.range(0);
        auto mat = DenseMatrix<T>::template random_uniform<RandomSource>(size, size);
        const auto vec = VectorND<T>::template random_uniform<RandomSource>(size);
        for (auto _ : state) {
            if constexpr (Left)
                PHYSICA_BENCH(applyHouseholder(T(2), vec, mat));
            else
                PHYSICA_BENCH(applyHouseholder(mat, T(2), vec));
            benchmark::DoNotOptimize(mat);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(kernel<true>)->Name("applyHouseHolder left")->Arg(512);
BENCHMARK(kernel<false>)->Name("applyHouseHolder right")->Arg(512);
