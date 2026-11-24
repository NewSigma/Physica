/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Scalar/Diff.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void forward(float64 v, float64& g) noexcept {
        using dfloat = Diff<float64, DiffMode::Forward>;
        dfloat x(v, 1);
        dfloat y = sin(x) / x;
        g = y.grad();
    }

    void reverse(float64 v, float64& g) noexcept {
        using dfloat = Diff<float64, DiffMode::Reverse>;
        dfloat x(v);
        (sin(x) / x).reverse();
        g = x.grad();
    }
}

BENCHMARK([](benchmark::State& state) {
    RandomSource::getInstance().reseed();

    float64 v = float64::random_uniform<RandomSource>();
    float64 g;
    for (auto _ : state) {
        PHYSICA_BENCH(forward(v, g));
        benchmark::DoNotOptimize(g);
    }
})->Name("Diff forward");
BENCHMARK([](benchmark::State& state) {
    RandomSource::getInstance().reseed();

    float64 v = float64::random_uniform<RandomSource>();
    float64 g;
    for (auto _ : state) {
        PHYSICA_BENCH(reverse(v, g));
        benchmark::DoNotOptimize(g);
    }
})->Name("Diff reverse");
