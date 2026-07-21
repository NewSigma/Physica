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
#include <cmath>
#include <chrono>
#include <fstream>
#include <print>
#include <unordered_map>
#include "Physica/Config.h"
#include "Benchmark.h"

import Physica.Core;

using namespace benchmark;
using namespace Physica;

namespace {
    constexpr std::array<std::size_t, 6> CacheSizes{
        HostDevAttr::LineSizeL1D / 2,
        HostDevAttr::LineSizeL1D,
        HostDevAttr::CacheSizeL1D * 3 / 4,
        HostDevAttr::CacheSizeL2 * 3 / 4,
        HostDevAttr::CacheSizeL3 * 3 / 4,
        HostDevAttr::CacheSizeL3  * 6 / 5
    };

    class Reporter : public ConsoleReporter {
        using Base = ConsoleReporter;
    public:
        Reporter();
        ~Reporter() = default;
        /* Operations */
        void PrintRunData(const Run& run) final;
    private:
        [[nodiscard]] static TimeUnit selectTimeUnit(const Run& r);
    };

    Reporter::Reporter() : ConsoleReporter(OutputOptions::OO_Tabular) {}

    void Reporter::PrintRunData(const Run& run) {
        auto r = run;
        r.time_unit = selectTimeUnit(r);
        Base::PrintRunData(r);
    }

    TimeUnit Reporter::selectTimeUnit(const Run& r) {
        const auto time = r.GetAdjustedRealTime();
        if (time > 1E9)
            return kSecond;
        if (time > 1E6)
            return kMillisecond;
        return time > 1E3 ? kMicrosecond : kNanosecond;
    }
}

namespace Physica {
    size_t makeVectorSize(int64_t level, size_t sizeElem) noexcept {
        return CacheSizes.at(level) / sizeElem;
    }

    size_t makeMatrixSize(int64_t level, size_t sizeElem) noexcept {
        return std::lround(std::sqrt(double(CacheSizes.at(level)) / double(sizeElem)));
    }
}

int main(int argc, char** argv) {
    if (!argv) {
        const char* const Executable = "Benchmark";
        argc = 1;
        argv = const_cast<char**>(&Executable);
    }
    Initialize(&argc, argv);
    if (ReportUnrecognizedArguments(argc, argv))
        return 1;

    std::ofstream("Version", std::ios_base::trunc) << Physica::version();
    auto start = std::chrono::steady_clock::now();
    {
        Reporter reporter{};
        RunSpecifiedBenchmarks(&reporter);
        Shutdown();
    }
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start);
    std::println("Elapsed time: {} Min", duration.count());
    return 0;
}
