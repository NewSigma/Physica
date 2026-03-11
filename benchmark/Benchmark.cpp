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
#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <print>
#include <unordered_map>
#include "Physica/Core/Version.h"
#include "Physica/Config.h"
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
        HostDevAttr::CacheSizeL3  * 6 / 5
    };

    class Reporter : public ConsoleReporter {
        using Base = ConsoleReporter;
    private:
        std::unordered_map<std::string, std::vector<Run>> subReportMap;
    public:
        Reporter();
        ~Reporter() = default;
        /* Operations */
        bool ReportContext(const Context&) override { return true; }
        void ReportRuns(const std::vector<Run>& reports) override;
        void Finalize() override;
    private:
        [[nodiscard]] static TimeUnit selectTimeUnit(const Run& r);
    };

    Reporter::Reporter() : ConsoleReporter(OutputOptions::OO_None) {
        std::ofstream fout("Version", std::ios_base::trunc);
        fout << Physica::version();
    }

    void Reporter::ReportRuns(const std::vector<Run>& reports) {
        for (const auto& r : reports) {
            const std::string_view name = r.run_name.function_name;
            const size_t pos = name.find_first_not_of(' ');
            const size_t n = name.find_first_of(' ', pos) - pos;
            const auto first_word = std::string(name.substr(pos, n));

            subReportMap[first_word].push_back(r);
        }
    }

    void Reporter::Finalize() {
        for (auto& pair : subReportMap) {
            std::ofstream fout(pair.first, std::ios_base::trunc);
            Base::SetOutputStream(&fout);

            auto& subReports = pair.second;
            Base::PrintHeader(subReports[0]);
            for (auto& r : subReports) {
                r.time_unit = selectTimeUnit(r);
                Base::PrintRunData(r);
            }

        }
        Base::SetOutputStream(&std::cout);
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
    std::string makeBenchID(std::source_location loc) noexcept {
        std::string result = loc.file_name();
        {
            std::string cur = std::source_location::current().file_name();
            size_t beginBaseDir = cur.find_last_of('/');
            result = result.substr(beginBaseDir + 1);
        }
        {
            size_t lastDot = result.find_last_of('.');
            assert(lastDot != std::string::npos && "[Error]: Unexpected missing extension");
            result = result.substr(0, lastDot);
        }
        std::ranges::replace(result, '/', '.');
        return result;
    }

    size_t makeVectorSize(int64_t level, size_t sizeElem) noexcept {
        return CacheSizes.at(level) / sizeElem;
    }

    size_t makeMatrixSize(int64_t level, size_t sizeElem) noexcept {
        return std::lround(std::sqrt(double(CacheSizes.at(level)) / double(sizeElem)));
    }
}

int main(int argc, char** argv) {
    const char* const Executable = "Benchmark";
    if (!argv) {
        argc = 1;
        argv = const_cast<char**>(&Executable);
    }
    Initialize(&argc, argv);
    if (ReportUnrecognizedArguments(argc, argv))
        return 1;

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
