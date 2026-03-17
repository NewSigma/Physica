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
#include "Physica/Core/Scalar/Diff.h"
#include "Benchmark.h"

using namespace Physica;

namespace {
    void forward(float64& g) noexcept {
        using dfloat = Diff<float64, DiffMode::Forward>;
        dfloat x(0.5, 1);
        dfloat y = sin(x) / x;
        g = y.grad();
    }

    void reverse(float64& g) noexcept {
        using dfloat = Diff<float64, DiffMode::Reverse>;
        dfloat x(0.5);
        (sin(x) / x).reverse();
        g = x.grad();
    }
}

BENCHMARK([](benchmark::State& state) {
    float64 g;
    for (auto _ : state) {
        PHYSICA_BENCH(forward(g));
        benchmark::DoNotOptimize(g);
    }
})->Name("Diff forward");
BENCHMARK([](benchmark::State& state) {
    float64 g;
    for (auto _ : state) {
        PHYSICA_BENCH(reverse(g));
        benchmark::DoNotOptimize(g);
    }
})->Name("Diff reverse");
