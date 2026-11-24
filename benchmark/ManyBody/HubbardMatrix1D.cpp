/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/KFermiRepr.h"
#include "Benchmark.h"

using namespace Physica;
using RealType = float64;
using ScalarType = Complex<RealType>;
using RandomSource = Random<MCG, std::mt19937::default_seed>;

namespace {
    void func(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        constexpr unsigned int NumSite = 10;
        constexpr unsigned int NumParticle = NumSite / 2;
        constexpr double HoppingT = 1.0;
        constexpr double RepelU = 4;
        using ReprType = KFermiRepr<1, NumSite, true>;

        const SquareLattice<1> lattice({NumSite}, 1);
        const HubbardMatrix<ScalarType, ReprType> model(HoppingT, RepelU, lattice, ReprType({NumParticle, NumParticle}, 0));
        auto v = VectorND<ScalarType>::random_uniform<RandomSource>(model.getRow());
        auto expr = model * v;

        VectorND<ScalarType> v1(model.getRow());
        for (auto _ : state) {
            PHYSICA_BENCH(expr.assign(v1));
            benchmark::DoNotOptimize(v1);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(func)->Name("HubbardMatrix1D");
