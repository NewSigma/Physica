/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <cmath>
#include <fstream>
#include <print>
#include "Physica/Core/Version.h"
#include "Physica/Config.h"
#include "Utils/Timer.h"
#include "Benchmark.h"

using namespace benchmark;
using namespace Physica;

namespace {
    constexpr std::array<std::size_t, 6> CacheSizes{
        HostDevAttr::LineSizeL1D / 2,
        HostDevAttr::LineSizeL1D,
        HostDevAttr::CacheSizeL1D * 3 / 4,
        HostDevAttr::CacheSizeL2 * 3 / 4,
        HostDevAttr::CacheSizeL3 * 3 / 4,
        HostDevAttr::CacheSizeL3 * 6 / 5,
    };
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
    HostTimer timer;
    timer.tic();
    RunSpecifiedBenchmarks();
    Shutdown();
    std::println(stderr, "Elapsed time: {} Min", static_cast<int64_t>(timer.toc() / 60.0));
    return 0;
}
