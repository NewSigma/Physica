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
#include <gperftools/profiler.h>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"
#include "Physica/Core/Physics/ManyBody/TPQ.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomSource = Random<>;
constexpr unsigned int NumSiteX = 4;
constexpr unsigned int NumSiteY = 2;
constexpr unsigned int NumSite = NumSiteX * NumSiteY;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr double Beta = 16;

namespace {
    void bench(benchmark::State& state) {
        using ReprType = FermiRepr<2, NumSite, true>;
        using Hamilton = HubbardMatrix<ScalarType, ReprType>;
        const SquareLattice<2> lattice({NumSiteX, NumSiteY}, 1);
        const Hamilton hamilton(HoppingT, RepelU, lattice, ReprType(4, 4));
        auto psi = TPQ<ScalarType>::random_normal<RandomSource>(hamilton.getNumState());
        psi.pre_nvt_step(hamilton, Beta);
        psi.random_normal<RandomSource>();
        for (auto _ : state)
            [[clang::noinline]] psi.nvt_step(hamilton, Beta);
    }
}

BENCHMARK(bench)->Name("TPQ");
