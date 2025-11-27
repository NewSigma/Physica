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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"

using namespace Physica;

namespace {
    template<bool NeedEigenVec>
    void direct(benchmark::State& state) {
        using T = float64;
        const auto size = state.range(0);
        SymmEigenSolver<T> solver(size, NeedEigenVec);

        auto m = MatrixND<T>::random_uniform<Random<MT19937>>(size);
        for (auto _ : state) {
            solver.compute(m);
            benchmark::DoNotOptimize(solver);
            benchmark::ClobberMemory();
        }
    }

    template<bool NeedEigenVec>
    void base(benchmark::State& state) {
        using T = float64;
        const auto size = state.range(0);
        SymmEigenSolver<T> solver(size, NeedEigenVec);

        auto m = MatrixND<T>::random_uniform<Random<MT19937>>(size);
        for (auto _ : state) {
            solver.compute_base(m);
            benchmark::DoNotOptimize(solver);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(direct<false>)->Name("SymmEigenSolver s")->Arg(4)->Arg(64)->Arg(1024);
BENCHMARK(direct<true>)->Name("SymmEigenSolver sv")->Arg(4)->Arg(64)->Arg(1024);
BENCHMARK(base<false>)->Name("SymmEigenSolver s base")->Arg(4)->Arg(64)->Arg(1024);
BENCHMARK(base<true>)->Name("SymmEigenSolver sv base")->Arg(4)->Arg(64)->Arg(1024);
