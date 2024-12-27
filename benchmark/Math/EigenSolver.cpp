/*
 * Copyright 2024 Weibo He.
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

template<bool EigVector>
static void direct(benchmark::State& state) {
    using T = float64;
    using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element>;
    const MatrixType data = TransIsingMatrix(TransIsing<T, 1>({{10}, 1}, 1, 0.01), SpinRepr<1, 10>(10));
    EigenSolver<T> solver(data.getRow());
    for (auto _ : state)
        solver.compute(data, EigVector);
}

template<bool EigVector>
static void base(benchmark::State& state) {
    using T = float64;
    using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element>;
    const MatrixType data = TransIsingMatrix(TransIsing<T, 1>({{10}, 1}, 1, 0.01), SpinRepr<1, 10>(10));
    EigenSolver<T> solver(data.getRow());
    for (auto _ : state)
        solver.compute_base(data, EigVector);
}

BENCHMARK(direct<false>)->Name("EigenSolver_s_direct")->Unit(benchmark::kMillisecond);
BENCHMARK(direct<true>)->Name("EigenSolver_sv_direct")->Unit(benchmark::kMillisecond);
BENCHMARK(base<false>)->Name("EigenSolver_s_base")->Unit(benchmark::kMillisecond);
BENCHMARK(base<true>)->Name("EigenSolver_sv_base")->Unit(benchmark::kMillisecond);
