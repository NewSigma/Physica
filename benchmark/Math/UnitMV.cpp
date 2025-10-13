/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void kernel1(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        const auto v1 = VectorND<T>::random_uniform<RandomSource>(size);
        const auto v2 = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            buffer = v1 - v2 * T(2.0);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void kernel2(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        const auto v1 = VectorND<T>::random_uniform<RandomSource>(size);
        const auto v2 = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            buffer = v1 - (UnitMatrix<T>(size) * T(2.0)) * v2;
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void kernel3(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        const auto v = VectorND<T>::random_uniform<RandomSource>(size);
        VectorND<T> buffer(size);
        for (auto _ : state) {
            buffer += (UnitMatrix<T>(size) * T(2.0)) * v;
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(kernel1)->Name("UnitMV 1")->Arg(256);
BENCHMARK(kernel2)->Name("UnitMV 2")->Arg(256);
BENCHMARK(kernel3)->Name("UnitMV 3")->Arg(256);
