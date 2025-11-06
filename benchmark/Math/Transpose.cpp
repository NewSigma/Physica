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
#include <benchmark/benchmark.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void trans(benchmark::State& state) {
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixOption::Col>;
        const size_t order = state.range(0);
        const auto m = MatrixType::template random_uniform<RandomSource>(order, order);
        MatrixType buffer(order, order);
        for (auto _ : state) {
            m.transpose().assign(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void trans_base(benchmark::State& state) {
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixOption::Col>;
        const size_t order = state.range(0);
        const auto m = MatrixType::template random_uniform<RandomSource>(order, order);
        MatrixType buffer(order, order);
        for (auto _ : state) {
            m.transpose().assign_base(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(trans)->Name("trans")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64);

BENCHMARK(trans_base)->Name("trans base")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64);
