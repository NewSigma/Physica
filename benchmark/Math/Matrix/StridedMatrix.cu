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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/CUDAEvent.cuh"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<>;
using T = float32;

namespace {
    enum class Region : int8_t { Full, Quarter };

    template<Matrix M, Region R>
    void assign(benchmark::State& state) {
        const size_t n = state.range(0);
        const size_t len = R == Region::Quarter ? n / 2 : n;
        auto src = device_obj<M>::template random_uniform<RandomSource>(n, n);
        auto dst = device_obj<M>::template random_uniform<RandomSource>(n, n);
        auto s = src.block(0, len, 0, len);
        auto d = dst.block(0, len, 0, len);
        CUDAExecutor::wait();

        auto tic = CUDAEvent();
        auto toc = CUDAEvent();
        for (auto _ : state) {
            tic.record();
            s.assign(d);
            toc.record();
            toc.wait();

            float elapsedTime{};
            cudaEventElapsedTime(&elapsedTime, tic, toc);
            state.SetIterationTime(elapsedTime / 1000);
        }
    }
}

using ColMatrix = DenseMatrix<T, MatrixMajor::Col>;
using RowMatrix = DenseMatrix<T, MatrixMajor::Row>;

BENCHMARK_TEMPLATE(assign, RowMatrix, Region::Full)
        ->Name("strided assign row full")->UseManualTime()->Arg(64)->Arg(256)->Arg(1024)->Arg(4096);
BENCHMARK_TEMPLATE(assign, RowMatrix, Region::Quarter)
        ->Name("strided assign row quarter")->UseManualTime()->Arg(64)->Arg(256)->Arg(1024)->Arg(4096);
BENCHMARK_TEMPLATE(assign, ColMatrix, Region::Full)
        ->Name("strided assign col full")->UseManualTime()->Arg(64)->Arg(256)->Arg(1024)->Arg(4096);
BENCHMARK_TEMPLATE(assign, ColMatrix, Region::Quarter)
        ->Name("strided assign col quarter")->UseManualTime()->Arg(64)->Arg(256)->Arg(1024)->Arg(4096);
