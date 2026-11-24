/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Utils/Timer.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<PCG64DXSM, std::mt19937::default_seed>;
using T = float32;

namespace {
    enum class Region : int8_t { Full, Quarter };

    template<Matrix M, Region R>
    void assign(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t n = state.range(0);
        const size_t len = R == Region::Quarter ? n / 2 : n;
        auto src = device_obj<M>::template random_uniform<RandomSource>(n, n);
        auto dst = device_obj<M>::template random_uniform<RandomSource>(n, n);
        auto s = src.block(0, len, 0, len);
        auto d = dst.block(0, len, 0, len);
        CUDAExecutor::wait();

        DeviceTimer timer;
        for (auto _ : state) {
            timer.tic();
            s.assign(d);
            state.SetIterationTime(timer.toc());
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
