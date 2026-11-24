/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/TransIsingMatrix.h"
#include "Benchmark.h"

using namespace Physica;

namespace {
    template<bool NeedEigenVec>
    void base(benchmark::State& state) {
        using T = float64;
        const SquareLattice<1> lattice({10}, 1);
        const DenseMatrix<T> data = TransIsingMatrix<T, 1, 10>(1, 0.01, lattice);
        EigenSolver<T> solver(data.getRow(), NeedEigenVec);
        for (auto _ : state) {
            PHYSICA_BENCH(solver.compute_base(data));
            benchmark::DoNotOptimize(solver);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<bool NeedEigenVec>
    void base_compare(benchmark::State& state) {
        using T = float64;
        const SquareLattice<1> lattice({10}, 1);
        const DenseMatrix<T> data = TransIsingMatrix<T, 1, 10>(1, 0.01, lattice);
        EigenSolver<T> solver(data.getRow(), NeedEigenVec);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> EigenSolver<T>& { solver.compute_base(data); return solver; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> EigenSolver<T>& { solver.compute_mkl(data); return solver; }, "mkl"),
                        backend([&](auto&&...) -> EigenSolver<T>& { solver.compute(data); return solver; }, "dispatch"));
    }
#endif
}

BENCHMARK(base<false>)->Name("EigenSolver s base");
BENCHMARK(base<true>)->Name("EigenSolver sv base");

#ifndef PHYSICA_EmitLLVM
BENCHMARK(base_compare<false>)->Name("EigenSolver s")->UseManualTime();
BENCHMARK(base_compare<true>)->Name("EigenSolver sv")->UseManualTime();
#endif
