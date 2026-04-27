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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/TransIsingMatrix.h"
#include "Benchmark.h"

using namespace Physica;

namespace {
    template<bool NeedEigenVec>
    void direct(benchmark::State& state) {
        using T = float64;
        const SquareLattice<1> lattice({10}, 1);
        const DenseMatrix<T> data = TransIsingMatrix<T, 1, 10>(1, 0.01, lattice);
        EigenSolver<T> solver(data.getRow(), NeedEigenVec);
        for (auto _ : state) {
            solver.compute(data);
            benchmark::DoNotOptimize(solver);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(direct<false>)->Name("EigenSolver s default");
BENCHMARK(direct<true>)->Name("EigenSolver sv default");
