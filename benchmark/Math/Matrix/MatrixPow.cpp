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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixPow.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    void kernel(benchmark::State& state) {
        using T = float64;
        const auto size = state.range(0);
        const auto m = MatrixND<T>::random_uniform<RandomSource>(size);
        const auto v = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = pow(m, 3) * v;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
        }
    }
}

BENCHMARK(kernel)->Name("MatrixPow")->Arg(1024);
