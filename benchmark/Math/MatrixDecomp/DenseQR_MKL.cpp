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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseQR.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    void qr_mkl(benchmark::State& state) {
        using T = float64;
        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        DenseQR<T> qr(order, order);
        for (auto _ : state) {
            qr.compute_mkl(m);
            benchmark::DoNotOptimize(qr);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(qr_mkl)->Name("DenseQR mkl")
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024);
