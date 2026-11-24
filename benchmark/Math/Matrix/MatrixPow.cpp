/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixPow.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void kernel(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const auto size = state.range(0);
        const auto m = MatrixND<T>::random_uniform<RandomSource>(size);
        const auto v = VectorND<T>::random_uniform<RandomSource>(size);
        auto expr = pow(m, 3) * v;
        VectorND<T> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
        }
    }
}

BENCHMARK(kernel)->Name("MatrixPow")->Arg(1024);
