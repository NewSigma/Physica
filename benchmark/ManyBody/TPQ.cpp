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
#include <iostream>
#include <benchmark/benchmark.h>
#include <gperftools/profiler.h>
#include <Physica/Core/Math/Random/Random.h>
#include <Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h>
#include <Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h>
#include <Physica/Core/Physics/ManyBody/TPQ.h>
#include <Physica/Core/Utils/BenchmarkHelper.h>

using namespace Physica::Core;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomType = Random<std::mt19937>;
constexpr unsigned int NumSiteX = 4;
constexpr unsigned int NumSiteY = 2;
constexpr unsigned int NumSite = NumSiteX * NumSiteY;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr double Beta = 16;

namespace {
    static void main(benchmark::State& state) {
        using ReprType = SpinRepr<2, NumSite, true>;
        using Hamilton = HubbardMatrix<ScalarType, ReprType>;
        const LatticeModel<2> lattice({NumSiteX, NumSiteY}, 1);
        const Hubbard<ScalarType, 2> hubbard(lattice, HoppingT, RepelU);
        const Hamilton hamilton(hubbard, ReprType(4, 4));
        auto& gen = RandomType::getInstance();
        auto psi = TPQ<ScalarType>::random_normal(hamilton.getNumState(), gen);
        psi.pre_nvt_step(hamilton, Beta);
        for (auto _ : state) {
            psi.random_normal(gen);
            psi.template nvt_step<Hamilton, SequentialExecutor>(hamilton, Beta);
            [[maybe_unused]] auto xi = psi.lnPartitionXi();
        };
    }
}

BENCHMARK(main)->Name("TPQ")->Unit(benchmark::kSecond);
