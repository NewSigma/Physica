/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Benchmark.h"

using namespace Physica;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void assign(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using T = float64;
        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        MatrixND<T> buffer(order);
        for (auto _ : state) {
            m.tril_unit().assign_base(buffer);
            benchmark::DoNotOptimize(buffer);
            benchmark::ClobberMemory();
        }
    }

#ifndef PHYSICA_EmitLLVM
    template<Scalar T>
    void assign_compare(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        const size_t order = state.range(0);
        const auto m = MatrixND<T>::template random_uniform<RandomSource>(order, order);
        auto expr = m.tril_unit();
        MatrixND<T> buffer(order);
        rotate_backends(state, HostTimer{},
                        backend([&](auto&&...) -> MatrixND<T>& { expr.assign_base(buffer); return buffer; }, "base"),
                        backend<HasMKL()>([&](auto&&...) -> MatrixND<T>& { expr.assign_mkl(buffer); return buffer; }, "mkl"),
                        backend([&](auto&&...) -> MatrixND<T>& { expr.assign(buffer); return buffer; }, "dispatch"));
    }
#endif
}

BENCHMARK(assign)->Name("Trig assign base")
    ->Arg(256)
    ->Arg(1024)
    ->Arg(8192);

#ifndef PHYSICA_EmitLLVM
BENCHMARK(assign_compare<float64>)->Name("Trig assign")->UseManualTime()->Arg(256)->Arg(1024)->Arg(8192);
#endif
