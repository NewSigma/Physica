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
#include "Physica/CRCoro.h"
#include "Benchmark.h"

namespace {
    struct Int : public Physica::CRCoro<Int> {
        int64_t i;

        Int(int64_t i_) : i(i_) {}
    };

    int64_t sum1(int64_t n) noexcept {
        if (n == 0)
            return 0;
        return n + sum1(n - 1);
    }

    Int sum2(int64_t n) noexcept {
        if (n == 0)
            co_return 0;
        co_return n + sum2(n - 1).i;
    }

    void bench1(benchmark::State &state) {
        int64_t n = state.range(0);
        int64_t a{};
        for (auto _ : state) {
            PHYSICA_BENCH(a = sum1(n));
            benchmark::DoNotOptimize(a);
        }
    }

    void bench2(benchmark::State &state) {
        int64_t n = state.range(0);
        int64_t a{};
        for (auto _ : state) {
            PHYSICA_BENCH(a = sum2(n).i);
            benchmark::DoNotOptimize(a);
        }
    }
}

BENCHMARK(bench1)->Name("CRCoro base")->Arg(2048);
BENCHMARK(bench2)->Name("CRCoro")->Arg(2048);