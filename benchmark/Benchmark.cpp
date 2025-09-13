/*
 * Copyright 2024-2025 Weibo He.
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
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <benchmark/benchmark.h>
#include "Physica/Core/Version.h"

using namespace benchmark;

namespace {
    const char* const Executable = "Benchmark";

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
        for (const auto& pair : subReportMap) {
            std::ofstream fout(pair.first, std::ios_base::trunc);
            Base::SetOutputStream(&fout);

            const auto& subReports = pair.second;
            Base::PrintHeader(subReports[0]);
            for (const auto& r : subReports)
                Base::PrintRunData(r);
        }
        Base::SetOutputStream(&std::cout);
    }
}

int main(int argc, char** argv) {
    if (!argv) {
        argc = 1;
        argv = const_cast<char**>(&Executable);
    }
    Initialize(&argc, argv);
    if (ReportUnrecognizedArguments(argc, argv))
        return 1;

    Reporter reporter{};
    RunSpecifiedBenchmarks(&reporter);
    Shutdown();
    return 0;
}
