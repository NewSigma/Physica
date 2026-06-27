/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Parallel/ThreadPool.h"
#include "Physica/Core/Utils/Cycler.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    // scheduling latency
    void bench_schedule(benchmark::State& state) {
        // Warmup
        std::ignore = Cycler::getCyclesPerSec();
        auto& pool = ThreadPool::getInstance();
        pool.restart();

        auto x = float64::random_uniform<RandomSource>();
        for (auto _ : state) {
            auto from = Cycler::tic();
            PHYSICA_BENCH(auto task = schedule<Thread>([&x]() noexcept {
                x = reciprocal(x);
            }));
            auto to = Cycler::toc();
            state.SetIterationTime(Cycler::toSeconds(to - from));
            task.wait();
        }
        pool.shouldExit();
    }

    void bench_parallel_for1(benchmark::State& state) {
        // Warmup
        std::ignore = Cycler::getCyclesPerSec();
        auto& pool = ThreadPool::getInstance();
        pool.restart();

        auto x = VectorND<float64>::random_uniform<RandomSource>(pool.getNumThreads());
        for (auto _ : state) {
            auto from = Cycler::tic();
            PHYSICA_BENCH(auto task = parallel_for<Thread>([&x](size_t i) noexcept {
                x[i] = reciprocal(x[i]);
            }, x.getLength()));
            auto to = Cycler::toc();
            state.SetIterationTime(Cycler::toSeconds(to - from));
            task.wait();
        }
        pool.shouldExit();
    }

    void bench_parallel_for2(benchmark::State& state) {
        // Warmup
        std::ignore = Cycler::getCyclesPerSec();
        auto& pool = ThreadPool::getInstance();
        pool.restart();

        auto x = VectorND<float64>::random_uniform<RandomSource>(1024);
        for (auto _ : state) {
            auto from = Cycler::tic();
            PHYSICA_BENCH(auto task = parallel_for<Thread>([&x](size_t i) noexcept {
                x[i] = reciprocal(x[i]);
            }, x.getLength(), pool.getNumThreads()));
            auto to = Cycler::toc();
            state.SetIterationTime(Cycler::toSeconds(to - from));
            task.wait();
        }
        pool.shouldExit();
    }
}

BENCHMARK(bench_schedule)->Name("Thread schecule")->UseManualTime();
BENCHMARK(bench_parallel_for1)->Name("Thread parallel_for 1")->UseManualTime();
BENCHMARK(bench_parallel_for2)->Name("Thread parallel_for 2")->UseManualTime();
