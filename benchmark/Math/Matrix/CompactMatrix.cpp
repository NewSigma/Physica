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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<MCG>;

namespace {
    // Baseline for mismatch
    void assign_match(benchmark::State& state) {
        const auto size = makeMatrixSize(state.range(0), sizeof(T));
        const MatrixND<T> data = MatrixND<T>::template random_uniform<RandomSource>(size);
        MatrixND<T> buffer(size);
        for (auto _ : state) {
            data.assign(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void assign_mismatch(benchmark::State& state) {
        const auto size = makeMatrixSize(state.range(0), sizeof(T));
        const auto data = DenseMatrix<T, MatrixMajor::Col>::template random_uniform<RandomSource>(size);
        DenseMatrix<T, MatrixMajor::Row> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(data.assign_base(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(assign_match)->Name("Compact match")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(assign_mismatch)->Name("Compact mismatch base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
