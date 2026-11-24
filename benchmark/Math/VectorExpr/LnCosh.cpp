/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    template<Scalar T>
    void lncosh(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const VectorND<T> data = VectorND<T>::template random_uniform<RandomSource>(1024);
        auto expr = lncosh(data);
        VectorND<T> buffer(1024);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(lncosh<float32>)->Name("lncosh float32");
BENCHMARK(lncosh<float64>)->Name("lncosh float64");
BENCHMARK(lncosh<cfloat32>)->Name("lncosh cfloat32");
BENCHMARK(lncosh<cfloat64>)->Name("lncosh cfloat64");
