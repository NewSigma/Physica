/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <benchmark/benchmark.h>
#include "Benchmark.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Utils/Timer.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;
using T = float32;
constexpr int N = 1024;

namespace {
    template<int Major>
    void assign_base(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        auto A = device_obj<DenseMatrix<T, Major>>::template random_uniform<RandomSource>(N, N);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(N);
        device_obj<VectorND<T>> y(N);

        DeviceTimer timer;
        for (auto _ : state) {
            timer.tic();
            (A * x).assign_base(y);
            state.SetIterationTime(timer.toc());
        }
    }

    template<int Major>
    void assign_cublas(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        auto A = device_obj<DenseMatrix<T, Major>>::template random_uniform<RandomSource>(N, N);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(N);
        device_obj<VectorND<T>> y(N);

        DeviceTimer timer;
        for (auto _ : state) {
            timer.tic();
            (A * x).assign_cublas(y);
            state.SetIterationTime(timer.toc());
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<int Major>
    void assign_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        auto A = device_obj<DenseMatrix<T, Major>>::template random_uniform<RandomSource>(N, N);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(N);
        device_obj<VectorND<T>> y(N);

        rotate_backends(state, DeviceTimer{},
                        backend([&](auto&&...) { (A * x).assign_base(y); }, "base"),
                        backend([&](auto&&...) { (A * x).assign_cublas(y); }, "cublas"));
    }
#endif
}

BENCHMARK(assign_base<MatrixMajor::Row>)->Name("GEMV cuda assign row base")->UseManualTime();
BENCHMARK(assign_base<MatrixMajor::Col>)->Name("GEMV cuda assign col base")->UseManualTime();
BENCHMARK(assign_cublas<MatrixMajor::Row>)->Name("GEMV cuda assign row cublas")->UseManualTime();
BENCHMARK(assign_cublas<MatrixMajor::Col>)->Name("GEMV cuda assign col cublas")->UseManualTime();

#ifndef PHYSICA_EmitLLVM
BENCHMARK(assign_compare<MatrixMajor::Row>)->Name("GEMV cuda assign row")->UseManualTime();
BENCHMARK(assign_compare<MatrixMajor::Col>)->Name("GEMV cuda assign col")->UseManualTime();
#endif
