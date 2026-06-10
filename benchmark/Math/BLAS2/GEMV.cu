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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Parallel/CUDAEvent.cuh"
#include "Physica/Core/Scalar/Complex.h"

using namespace Physica;
using RandomSource = Random<MCG>;
using T = float32;
constexpr int N = 1024;

namespace {
    template<int Major>
    void assign_base(benchmark::State& state) {
        auto A = device_obj<DenseMatrix<T, Major>>::template random_uniform<RandomSource>(N, N);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(N);
        device_obj<VectorND<T>> y(N);

        auto tic = CUDAEvent();
        auto toc = CUDAEvent();
        for (auto _ : state) {
            tic.record();
            (A * x).assign_base(y);
            toc.record();

            toc.wait();
            float elapsedTime{};
            cudaEventElapsedTime(&elapsedTime, tic, toc);
            state.SetIterationTime(elapsedTime / 1000);
        }
    }

    template<int Major>
    void assign_cublas(benchmark::State& state) {
        auto A = device_obj<DenseMatrix<T, Major>>::template random_uniform<RandomSource>(N, N);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(N);
        device_obj<VectorND<T>> y(N);

        auto tic = CUDAEvent();
        auto toc = CUDAEvent();
        for (auto _ : state) {
            tic.record();
            (A * x).assign_cublas(y);
            toc.record();

            toc.wait();
            float elapsedTime{};
            cudaEventElapsedTime(&elapsedTime, tic, toc);
            state.SetIterationTime(elapsedTime / 1000);
        }
    }
}

BENCHMARK(assign_base<MatrixMajor::Row>)->Name("GEMV cuda assign row base")->UseManualTime();
BENCHMARK(assign_base<MatrixMajor::Col>)->Name("GEMV cuda assign col base")->UseManualTime();
BENCHMARK(assign_cublas<MatrixMajor::Row>)->Name("GEMV cuda assign row cublas")->UseManualTime();
BENCHMARK(assign_cublas<MatrixMajor::Col>)->Name("GEMV cuda assign col cublas")->UseManualTime();
