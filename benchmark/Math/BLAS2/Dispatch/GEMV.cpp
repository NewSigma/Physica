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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    template<Scalar T, int Major>
    void gemv(benchmark::State& state) {
        const size_t size = state.range(0);
        const auto mat = DenseMatrix<T, Major>::template random_uniform<RandomSource>(size, size);
        const auto vec = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = mat * vec;
        VectorND<T> result(size);
        for (auto _ : state) {
            expr.assign(result);
            benchmark::DoNotOptimize(result);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(gemv<float64, MatrixMajor::Col>)->Name("GEMV col dispatch")->Arg(16)->Arg(32)->Arg(512);
BENCHMARK(gemv<float64, MatrixMajor::Row>)->Name("GEMV row dispatch")->Arg(16)->Arg(32)->Arg(512);
