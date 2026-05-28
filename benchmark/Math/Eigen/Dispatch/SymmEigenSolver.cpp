/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    template<bool NeedEigenVec>
    void direct(benchmark::State& state) {
        using T = float64;
        const auto size = state.range(0);
        SymmEigenSolver<T> solver(size, NeedEigenVec);

        auto m = MatrixND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            solver.compute(m);
            benchmark::DoNotOptimize(solver);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(direct<false>)->Name("SymmEigenSolver s dispatch")->Arg(4)->Arg(64)->Arg(1024);
BENCHMARK(direct<true>)->Name("SymmEigenSolver sv dispatch")->Arg(4)->Arg(64)->Arg(1024);
