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
#include <benchmark/benchmark.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/CUDAEvent.cuh"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"

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

BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::Full)
        ->Name("strided assign full")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::Minor)
        ->Name("strided assign minor")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::Mid)
        ->Name("strided assign mid")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK_TEMPLATE(assign, Tensor3D<T>, Slice::OneEighth)
        ->Name("strided assign one_eighth")->UseManualTime()->Arg(32)->Arg(64)->Arg(128)->Arg(256);
