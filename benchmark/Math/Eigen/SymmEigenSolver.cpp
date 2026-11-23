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
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    template<bool NeedEigenVec>
    void base(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const auto size = state.range(0);
        SymmEigenSolver<T> solver(size, NeedEigenVec);

        auto m = MatrixND<T>::random_uniform<RandomSource>(size);
        for (auto _ : state) {
            PHYSICA_BENCH(solver.compute_base(m));
            benchmark::DoNotOptimize(solver);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<bool NeedEigenVec>
    void base_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const auto size = state.range(0);
        SymmEigenSolver<T> solver(size, NeedEigenVec);
        auto m = MatrixND<T>::random_uniform<RandomSource>(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> SymmEigenSolver<T>& { solver.compute_base(m); return solver; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> SymmEigenSolver<T>& { solver.compute_mkl(m); return solver; }, "mkl"),
                        backend([&](auto&&...) -> SymmEigenSolver<T>& { solver.compute(m); return solver; }, "dispatch"));
    }
#endif
}

BENCHMARK(base<false>)->Name("SymmEigenSolver s base")->Arg(4)->Arg(64)->Arg(1024);
BENCHMARK(base<true>)->Name("SymmEigenSolver sv base")->Arg(4)->Arg(64)->Arg(1024);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(base_compare<false>)->Name("SymmEigenSolver s")->UseManualTime()->Arg(4)->Arg(64)->Arg(1024);
BENCHMARK(base_compare<true>)->Name("SymmEigenSolver sv")->UseManualTime()->Arg(4)->Arg(64)->Arg(1024);
#endif
