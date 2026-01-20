/*
 * Copyright 2024-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    void kernel1(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        const auto v1 = VectorND<T>::random_uniform<RandomSource>(size);
        const auto v2 = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = v1 - v2 * T(2.1);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            expr.assign(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void kernel2(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        const auto v1 = VectorND<T>::random_uniform<RandomSource>(size);
        const auto v2 = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = v1 - (IdentityMatrix<T>(size) * T(2.1)) * v2;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void kernel3(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        const auto v = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = (IdentityMatrix<T>(size) * T(2.1)) * v;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            expr.assign_add(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(kernel1)->Name("IdentityMV 1")->Arg(256);
BENCHMARK(kernel2)->Name("IdentityMV 2")->Arg(256);
BENCHMARK(kernel3)->Name("IdentityMV 3")->Arg(256);
