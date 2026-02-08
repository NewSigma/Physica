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
    void bench(benchmark::State& state) {
        // Warmup
        ThreadPool::numThreadRequired = Dynamic;
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

BENCHMARK(bench)->Name("ThreadPool")->UseManualTime();
