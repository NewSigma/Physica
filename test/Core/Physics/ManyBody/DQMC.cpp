/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/DQMC.h"
#include "Physica/Core/Physics/ManyBody/GreenSampler/ScalarSampler.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;
using RandomSource = Random<>;
constexpr int Dim = 2;

namespace {
    /**
     * Fuzzing test that half filling is free of sign problem
     */
    void halfFillTest() {
        constexpr double HoppingT = 1;
        constexpr double RepelU = 8;
        constexpr double Beta = 8;
        constexpr int NumSiteX = 4;
        constexpr int NumSiteY = 4;
        constexpr int NumSplit = Beta * 8;
        constexpr int NumSample = 1024;

        const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, NumSplit);
        auto dqmc = DQMC<T>(params);
        dqmc.step_random<RandomSource>();
        for (int i = 0; i < NumSample; ++i) {
            dqmc.step<RandomSource>();
            expect<RandomSource>(dqmc.getRSign().isPositive());
        }
    }
    // Test that trotter decomposition does not affect results of free system
    void freeFermion() {
        constexpr static T HoppingT = 1;
        constexpr static T RepelU = 0;
        constexpr static T Beta = 8;
        constexpr int NumSiteX = 8;
        constexpr int NumSiteY = 8;
        auto calcKinetic = [](int numSplit) static {
            const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
            const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, numSplit);
            auto dqmc = DQMC<T>(params);
            dqmc.step_random<RandomSource>();
            dqmc.step<RandomSource>();

            ScalarSampler<T> sampler(params, 1);
            sampler.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<T>::Kinetic);
            return sampler.calcMean();
        };

        constexpr double Prec = 1E-3;
        constexpr int N = 3;
        constexpr Array<int, N> splits{2, 4, 6};
        const auto kinetics = Array<T>::generate([&](size_t i) {
            return calcKinetic(splits[i]);
        }, N);

        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j)
                expect(scalarNear(kinetics[i], kinetics[j], Prec));
    }
}

int main() {
    halfFillTest();
    freeFermion();
    return 0;
}
