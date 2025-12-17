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
#include "Physica/Core/Physics/ManyBody/ReprSpace/KFermiRepr.h"

using namespace Physica;

namespace {
    void func(benchmark::State& state) {
        constexpr unsigned int NumSite = 10;
        constexpr unsigned int NumParticle = NumSite / 2;
        constexpr double HoppingT = 1.0;
        constexpr double RepelU = 4;

        using RealType = float64;
        using ScalarType = Complex<RealType>;
        using RandomSource = Random<MCG>;
        using ReprType = KFermiRepr<1, NumSite, true>;

        const SquareLattice<1> lattice({NumSite}, 1);
        const HubbardMatrix<ScalarType, ReprType> model(HoppingT, RepelU, lattice, ReprType({NumParticle, NumParticle}, 0));
        auto v = VectorND<ScalarType>::random_uniform<RandomSource>(model.getRow());
        auto expr = model * v;

        VectorND<ScalarType> v1(model.getRow());
        for (auto _ : state) {
            
            [[clang::noinline]] expr.assign(v1);
            benchmark::DoNotOptimize(v1);
            benchmark::ClobberMemory();
        }
    }
}

BENCHMARK(func)->Name("HubbardMatrix1D");
