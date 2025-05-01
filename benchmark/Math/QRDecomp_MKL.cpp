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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/QRDecomp.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

template<int Order = 0>
static void qr_mkl(benchmark::State& state) {
    static_assert(Order >= 0);
    using T = float64;
    using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Order, Order>;
    const size_t order = state.range(0);
    const auto m = MatrixType::template random_uniform<RandomSource>(order, order);
    QRDecomp<float64> qr(order, order);
    for (auto _ : state) {
        qr.compute_mkl(m);
        benchmark::DoNotOptimize(qr);
        benchmark::ClobberMemory();
    }
}

BENCHMARK(qr_mkl<2>)->Name("QR mkl")->Unit(benchmark::kNanosecond)->Arg(2);
BENCHMARK(qr_mkl<4>)->Name("QR mkl")->Unit(benchmark::kNanosecond)->Arg(4);
BENCHMARK(qr_mkl<8>)->Name("QR mkl")->Unit(benchmark::kNanosecond)->Arg(8);
BENCHMARK(qr_mkl<16>)->Name("QR mkl")->Unit(benchmark::kNanosecond)->Arg(16);
BENCHMARK(qr_mkl<32>)->Name("QR mkl")->Unit(benchmark::kNanosecond)->Arg(32);
BENCHMARK(qr_mkl)->Name("QR mkl")->Unit(benchmark::kMicrosecond)->Arg(64);
BENCHMARK(qr_mkl)->Name("QR mkl")->Unit(benchmark::kMillisecond)->Arg(256);
BENCHMARK(qr_mkl)->Name("QR mkl")->Unit(benchmark::kMillisecond)->Arg(1024);
