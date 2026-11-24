/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
