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
#include <benchmark/benchmark.h>
#include "Physica/Core/Version.h"

using namespace benchmark;

namespace {
    const char* const Executable = "Benchmark";

    class Reporter : public ConsoleReporter {
        using Base = ConsoleReporter;
    private:
        std::vector<Run> reports;
    public:
        Reporter() : ConsoleReporter(OutputOptions::OO_None) {}
        ~Reporter() = default;

        bool ReportContext(const Context&) override { return true; }
        void ReportRuns(const std::vector<Run>& group) override {
            reports.insert(reports.end(), group.begin(), group.end());
        }
        void Finalize() override {
            std::vector<Run> subReports{};
            std::string header{};
            for (const auto& r : reports) {
                const std::string_view name = r.run_name.function_name;
                const size_t pos = name.find_first_not_of(' ');
                const size_t n = name.find_first_of(' ', pos) - pos;
                const auto first_word = name.substr(pos, n);
                
                if (first_word == header) {
                    subReports.push_back(r);
                    continue;
                }
                
                ReportSubRuns(header, subReports);

                header = std::string(first_word);
                subReports.push_back(r);
            }
            ReportSubRuns(header, subReports);
        }
    private:
        void ReportSubRuns(const std::string& header, std::vector<Run>& subReports) {
            if (header.empty())
                return;

            std::ofstream fout(header, std::ios_base::trunc);
            Base::SetOutputStream(&fout);
            Base::ReportRuns(subReports);
            subReports.clear();
        }
    };
}

int main(int argc, char** argv) {
    if (!argv) {
        argc = 1;
        argv = const_cast<char**>(&Executable);
    }
    Initialize(&argc, argv);
    if (ReportUnrecognizedArguments(argc, argv))
        return 1;

    {
        std::ofstream fout("Version", std::ios_base::trunc);
        fout << Physica::version() << '\n';
    }
    Reporter reporter{};
    RunSpecifiedBenchmarks(&reporter);
    Shutdown();
    return 0;
}
