/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <benchmark/benchmark.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Utils/Timer.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;
using T = float32;

namespace {
    enum class Slice : int8_t {
        Full,
        Minor,
        Mid,
        OneEighth,
    };

    template<class TensorT, Slice S>
    void assign(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t n = state.range(0);
        const Index3D shape{n, n, n};
        Index3D count{};
        if constexpr (S == Slice::Full)
            count = shape;
        else if constexpr (S == Slice::Minor)
            count = Index3D{n, n, n / 2};
        else if constexpr (S == Slice::Mid)
            count = Index3D{n, n / 2, n};
        else
            count = Index3D{n / 2, n / 2, n / 2};
        auto src = device_obj<TensorT>::template random_uniform<RandomSource>(shape);
        auto dst = device_obj<TensorT>::template random_uniform<RandomSource>(shape);
        auto s = src.block(Index3D{}, count);
        auto d = dst.block(Index3D{}, count);
        CUDAExecutor::wait();

        DeviceTimer timer;
        for (auto _ : state) {
            timer.tic();
            s.assign(d);
            state.SetIterationTime(timer.toc());
        }
    }
}

BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::Full)
        ->Name("strided assign full")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::Minor)
        ->Name("strided assign minor")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::Mid)
        ->Name("strided assign mid")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::OneEighth)
        ->Name("strided assign one_eighth")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
