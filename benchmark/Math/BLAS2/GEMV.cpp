/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG>;

namespace {
    template<Scalar T, int Major>
    void gemv(benchmark::State& state) {
        const size_t size = state.range(0);
        const auto mat = DenseMatrix<T, Major>::template random_uniform<RandomSource>(size, size);
        const auto vec = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = mat * vec;
        VectorND<T> result(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_base(result));
            benchmark::DoNotOptimize(result);
            benchmark::ClobberMemory();
        }
    }

    void complex_real(benchmark::State& state) {
        using T = float64;
        using Tc = cfloat64;
        constexpr size_t size = 32;
        const auto mat = MatrixND<cfloat64>::template random_uniform<RandomSource>(size, size);
        const auto vec = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = mat * vec;
        VectorND<Tc> result(size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign_base(result));
            benchmark::DoNotOptimize(result);
            benchmark::ClobberMemory();
        }
    }

    // Benchmark we recursively lower GEMV::grad() without falling back to RValueVector::grad
    void forward2_grad(benchmark::State& state) {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Forward, 2>;
        constexpr int Size = 32;

        auto A = MatrixND<dfloat>(Size, Size);
        auto x = VectorND<dfloat>(Size);
        A.random_uniform<RandomSource>();
        x.random_uniform<RandomSource>();
        A.grads() = T(1);
        x.grads() = T(1);
        auto expr = square(A * x);
        VectorND<dfloat> buffer(Size);
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(gemv<float64, MatrixMajor::Col>)->Name("GEMV base")->Arg(16)->Arg(32)->Arg(512);
BENCHMARK(complex_real)->Name("GEMV complex-real");
BENCHMARK(forward2_grad)->Name("GEMV forward2");
