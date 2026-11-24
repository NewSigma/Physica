/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void trans(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixMajor::Col>;
        const size_t order = state.range(0);
        const auto m = MatrixType::template random_uniform<RandomSource>(order, order);
        auto trans = m.transpose();
        MatrixType buffer(order, order);
        for (auto _ : state) {
            PHYSICA_BENCH(trans.assign_base(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void trans_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using MatrixType = DenseMatrix<T, MatrixMajor::Col>;
        const size_t order = state.range(0);
        const auto m = MatrixType::template random_uniform<RandomSource>(order, order);
        auto expr = m.transpose();
        MatrixType buffer(order, order);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> MatrixType& { expr.assign_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> MatrixType& { expr.assign_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> MatrixType& { expr.assign(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(trans)->Name("trans base")
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(trans_compare<float64>)->Name("trans")->UseManualTime()->Arg(2)->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64);
#endif
