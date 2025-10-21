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
#include <benchmark/benchmark.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/TransIsingMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h"

using namespace Physica;

namespace {
    template<bool NeedEigenVec>
    void kernel(benchmark::State& state) {
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixOption::Col>;
        const SquareLattice<1> lattice({10}, 1);
        const MatrixType data = TransIsingMatrix<T, SpinRepr<1, 10>>(1, 0.01, lattice, SpinRepr<1, 10>(10));
        EigenSolver<T> solver(data.getRow(), NeedEigenVec);
        for (auto _ : state)
            solver.compute_mkl(data);
    }
}

BENCHMARK(kernel<false>)->Name("EigenSolver s mkl");
BENCHMARK(kernel<true>)->Name("EigenSolver sv mkl");
