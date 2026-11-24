/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"
#include "Physica/Core/Physics/ManyBody/TPQ.h"
#include "Benchmark.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomSource = Random<MCG, std::mt19937::default_seed>;
constexpr unsigned int NumSiteX = 4;
constexpr unsigned int NumSiteY = 2;
constexpr unsigned int NumSite = NumSiteX * NumSiteY;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr double Beta = 16;

namespace {
    void bench(benchmark::State& state) {
        RandomSource::getInstance().reseed();

        using ReprType = FermiRepr<2, NumSite, true>;
        using Hamilton = HubbardMatrix<ScalarType, ReprType>;
        const SquareLattice<2> lattice({NumSiteX, NumSiteY}, 1);
        const Hamilton hamilton(HoppingT, RepelU, lattice, ReprType(4, 4));
        auto psi = TPQ<ScalarType>::random_normal<RandomSource>(hamilton.getNumState());
        psi.pre_nvt_step(hamilton, Beta);
        psi.random_normal<RandomSource>();
        for (auto _ : state)
            PHYSICA_BENCH(psi.nvt_step(hamilton, Beta));
    }
}

BENCHMARK(bench)->Name("TPQ");
