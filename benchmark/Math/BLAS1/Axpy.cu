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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Physica/Core/Parallel/CUDAEvent.cuh"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    void axpy(benchmark::State& state) {
        using T = float32;
        const int64_t size = state.range(0);
        auto x = device_obj<VectorND<T>>::random_uniform<RandomSource>(size);
        device_obj<VectorND<T>> buffer(size, 0);
        auto expr = x * T(2);
        CUDAExecutor::wait();

        auto tic = CUDAEvent();
        auto toc = CUDAEvent();
        for (auto _ : state) {
            tic.record();
            expr.assign_add_base(buffer);
            toc.record();

            toc.wait();
            float elapsedTime{};
            cudaEventElapsedTime(&elapsedTime, tic, toc);
            state.SetIterationTime(elapsedTime / 1000);
        }
    }
}

constexpr int N = 1 << 20;
BENCHMARK(axpy)->Name("axpy cuda base")->UseManualTime()->Arg(N);
