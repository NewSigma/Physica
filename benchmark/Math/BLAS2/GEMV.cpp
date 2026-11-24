/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    template<Scalar T, int Major>
    void gemv(benchmark::State& state) {
        RandomSource::getInstance().reseed();

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
        RandomSource::getInstance().reseed();

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
        RandomSource::getInstance().reseed();

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

#ifndef PHYSICA_EmitLLVM
    template<Scalar T, int Major>
    void gemv_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t size = state.range(0);
        const auto mat = DenseMatrix<T, Major>::template random_uniform<RandomSource>(size, size);
        const auto vec = VectorND<T>::template random_uniform<RandomSource>(size);
        auto expr = mat * vec;
        VectorND<T> result(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign_base(result); return result; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> VectorND<T>& { expr.assign_mkl(result); return result; }, "mkl"),
                        backend([&](auto&&...) -> VectorND<T>& { expr.assign(result); return result; }, "dispatch"));
    }
#endif
}

BENCHMARK(gemv<float64, MatrixMajor::Col>)->Name("GEMV col base")->Arg(16)->Arg(32)->Arg(512);
BENCHMARK(gemv<float64, MatrixMajor::Row>)->Name("GEMV row base")->Arg(16)->Arg(32)->Arg(512);
BENCHMARK(complex_real)->Name("GEMV complex-real");
BENCHMARK(forward2_grad)->Name("GEMV forward2");

#ifndef PHYSICA_EmitLLVM
BENCHMARK(gemv_compare<float64, MatrixMajor::Col>)->Name("GEMV col")->UseManualTime()->Arg(16)->Arg(32)->Arg(512);
BENCHMARK(gemv_compare<float64, MatrixMajor::Row>)->Name("GEMV row")->UseManualTime()->Arg(16)->Arg(32)->Arg(512);
#endif
