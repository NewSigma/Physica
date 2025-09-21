/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

namespace {
    void scal(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        bool flag = true;
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            (x * T(flag ? std::numbers::pi : (1.0 / std::numbers::pi))).assign(x);
            benchmark::DoNotOptimize(x);
            benchmark::ClobberMemory();
            flag = !flag;
        }
    }

    void scal_base(benchmark::State& state) {
        using T = float64;
        const int64_t size = state.range(0);
        bool flag = true;
        auto x = VectorND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            (x * T(flag ? std::numbers::pi : (1.0 / std::numbers::pi))).assign_base(x);
            benchmark::DoNotOptimize(x);
            benchmark::ClobberMemory();
            flag = !flag;
        }
    }
}

BENCHMARK(scal)->Name("scal")->Arg(2);
BENCHMARK(scal)->Name("scal")->Arg(4);
BENCHMARK(scal)->Name("scal")->Arg(8);
BENCHMARK(scal)->Name("scal")->Arg(16);
BENCHMARK(scal)->Name("scal")->Arg(64);
BENCHMARK(scal)->Name("scal")->Arg(256);
BENCHMARK(scal)->Name("scal")->Arg(1024);
BENCHMARK(scal)->Name("scal")->Arg(16384);

BENCHMARK(scal_base)->Name("scal base")->Arg(2);
BENCHMARK(scal_base)->Name("scal base")->Arg(4);
BENCHMARK(scal_base)->Name("scal base")->Arg(8);
BENCHMARK(scal_base)->Name("scal base")->Arg(16);
BENCHMARK(scal_base)->Name("scal base")->Arg(64);
BENCHMARK(scal_base)->Name("scal base")->Arg(256);
BENCHMARK(scal_base)->Name("scal base")->Arg(1024);
BENCHMARK(scal_base)->Name("scal base")->Arg(16384);
