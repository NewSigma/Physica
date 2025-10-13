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
using RandomSource = Random<>;

namespace {
    template<int Order = 0>
    void qr(benchmark::State& state) {
        static_assert(Order >= 0);
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Order, Order>;
        const size_t order = state.range(0);
        const auto m = MatrixType::template random_uniform<RandomSource>(order, order);
        QRDecomp<float64> qr(order, order);
        for (auto _ : state) {
            qr.compute(m);
            benchmark::DoNotOptimize(qr);
            benchmark::ClobberMemory();
        }
    }

    template<int Order = 0>
    void qr_base(benchmark::State& state) {
        static_assert(Order >= 0);
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Order, Order>;
        const size_t order = state.range(0);
        const auto m = MatrixType::template random_uniform<RandomSource>(order, order);
        QRDecomp<float64> qr(order, order);
        for (auto _ : state) {
            qr.compute_base(m);
            benchmark::DoNotOptimize(qr);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(qr<2>)->Name("QR")->Arg(2);
BENCHMARK(qr<4>)->Name("QR")->Arg(4);
BENCHMARK(qr<8>)->Name("QR")->Arg(8);
BENCHMARK(qr<16>)->Name("QR")->Arg(16);
BENCHMARK(qr<32>)->Name("QR")->Arg(32);
BENCHMARK(qr)->Name("QR")
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024);

BENCHMARK(qr_base<2>)->Name("QR base")->Arg(2);
BENCHMARK(qr_base<4>)->Name("QR base")->Arg(4);
BENCHMARK(qr_base<8>)->Name("QR base")->Arg(8);
BENCHMARK(qr_base<16>)->Name("QR base")->Arg(16);
BENCHMARK(qr_base<32>)->Name("QR base")->Arg(32);
BENCHMARK(qr_base)->Name("QR base")
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024);
