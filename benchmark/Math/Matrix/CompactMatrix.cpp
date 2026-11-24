/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    // Baseline for mismatch
    void assign_match(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeMatrixSize(state.range(0), sizeof(T));
        const MatrixND<T> data = MatrixND<T>::template random_uniform<RandomSource>(size);
        MatrixND<T> buffer(size);
        for (auto _ : state) {
            data.assign(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

    void assign_mismatch(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeMatrixSize(state.range(0), sizeof(T));
        const auto data = DenseMatrix<T, MatrixMajor::Col>::template random_uniform<RandomSource>(size);
        DenseMatrix<T, MatrixMajor::Row> buffer(size);
        for (auto _ : state) {
            PHYSICA_BENCH(data.assign_base(buffer));
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void assign_mismatch_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const auto size = makeMatrixSize(state.range(0), sizeof(T));
        const auto data = DenseMatrix<T, MatrixMajor::Col>::template random_uniform<RandomSource>(size);
        DenseMatrix<T, MatrixMajor::Row> buffer(size);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> DenseMatrix<T, MatrixMajor::Row>& { data.assign_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> DenseMatrix<T, MatrixMajor::Row>& { data.assign_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> DenseMatrix<T, MatrixMajor::Row>& { data.assign(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(assign_match)->Name("Compact match")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
BENCHMARK(assign_mismatch)->Name("Compact mismatch base")->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(assign_mismatch_compare<float64>)->Name("Compact mismatch")->UseManualTime()->Arg(0)->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(5);
#endif
