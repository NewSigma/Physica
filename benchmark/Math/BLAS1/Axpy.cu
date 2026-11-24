/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Utils/Timer.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void axpy(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float32;
        const int64_t size = state.range(0);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(size);
        device_obj<VectorND<T>> buffer(size, 0);
        auto expr = x * T(2);
        CUDAExecutor::wait();

        DeviceTimer timer;
        for (auto _ : state) {
            timer.tic();
            expr.assign_add_base(buffer);
            state.SetIterationTime(timer.toc());
        }
    }

#ifndef PHYSICA_EmitLLVM
    void axpy_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float32;
        const int64_t size = state.range(0);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(size);
        device_obj<VectorND<T>> buffer(size, 0);
        auto expr = x * T(2);
        CUDAExecutor::wait();

        rotate_backends(state, DeviceTimer{},
                        backend([&](auto&&...) { expr.assign_add_base(buffer); }, "base"),
                        backend([&](auto&&...) { expr.assign_add_cublas(buffer); }, "cublas"),
                        backend([&](auto&&...) { expr.assign_add(buffer); }, "dispatch"));
    }
#endif
}

constexpr int N = 1 << 20;
BENCHMARK(axpy)->Name("axpy cuda base")->UseManualTime()->Arg(N);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(axpy_compare)->Name("axpy cuda")->UseManualTime()->Arg(N);
#endif
