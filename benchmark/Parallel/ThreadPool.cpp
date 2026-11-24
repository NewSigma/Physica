/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Parallel/Algorithm/Thread.h"
#include "Utils/Timer.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    // scheduling latency
    void bench_schedule(benchmark::State& state) {
        RandomSource::getInstance().reseed();
        std::ignore = ThreadPool::getInstance();

        HostTimer timer;
        auto x = float64::random_uniform<RandomSource>();
        for (auto _ : state) {
            timer.tic();
            PHYSICA_BENCH(auto task = schedule<Thread>([&x]() noexcept {
                x = reciprocal(x);
            }));
            state.SetIterationTime(timer.toc());
            task.wait();
        }
    }

    void bench_parallel_for1(benchmark::State& state) {
        RandomSource::getInstance().reseed();
        auto& pool = ThreadPool::getInstance();

        HostTimer timer;
        auto x = VectorND<float64>::random_uniform<RandomSource>(pool.getNumThreads());
        for (auto _ : state) {
            timer.tic();
            PHYSICA_BENCH(auto task = parallel_for<Thread>([&x](size_t i) noexcept {
                x[i] = reciprocal(x[i]);
            }, x.getLength()));
            state.SetIterationTime(timer.toc());
            task.wait();
        }
    }

    void bench_parallel_for2(benchmark::State& state) {
        RandomSource::getInstance().reseed();
        auto& pool = ThreadPool::getInstance();

        HostTimer timer;
        auto x = VectorND<float64>::random_uniform<RandomSource>(1024);
        for (auto _ : state) {
            timer.tic();
            PHYSICA_BENCH(auto task = parallel_for<Thread>([&x](size_t i) noexcept {
                x[i] = reciprocal(x[i]);
            }, x.getLength(), pool.getNumThreads()));
            state.SetIterationTime(timer.toc());
            task.wait();
        }
    }
}

BENCHMARK(bench_schedule)->Name("Thread schecule")->UseManualTime();
BENCHMARK(bench_parallel_for1)->Name("Thread parallel_for 1")->UseManualTime();
BENCHMARK(bench_parallel_for2)->Name("Thread parallel_for 2")->UseManualTime();
