/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/BayesDoping.h"
#include "Physica/Core/Physics/ManyBody/DQMC.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<PCG64DXSM, 5746352372098963138>;
constexpr int Dim = 2;
constexpr T HoppingT = 1;
constexpr T RepelU = 4;
constexpr T Beta = 4;
constexpr int NumSiteX = 4;
constexpr int NumSiteY = 4;
constexpr int NumSplit = Beta.toMachine() * 8;

namespace {
    void halfFilling() {
        // Test that half filling is trivially solved without effort
        constexpr int BayesIter = 4;
        const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
        HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, NumSplit);
        auto dqmc = DQMC<T>(params);
        auto dopant = BayesDoping<T>(Vegas<T, true>({T(-1), T(1)}, 10, 10), BayesIter, 1, 0, 4, 0);
        const T optimal = dopant.template solve<RandomSource>(T(1), params, dqmc);
        expect<RandomSource>(dopant.getNumSamples() == BayesIter + 1);
        expect<RandomSource>(optimal.isZero());
    }

    void doping() {
        constexpr T From = -4;
        constexpr T To = 0;
        constexpr int NumVegasSample = 10000;
        constexpr int NumBayesIter = 8;
        constexpr int NumVegasIter = 64;
        constexpr int NumWarmup = 32;
        constexpr int NumSample = 1024;
        constexpr int NumStepSGD = 100;
        constexpr T Target = 0.8;
        const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
        HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, NumSplit);
        auto dqmc = DQMC<T>(params);
        auto dopant = BayesDoping<T>(Vegas<T, true>({T(From), T(To)}, NumVegasSample, 30), NumBayesIter, NumVegasIter, NumWarmup, NumSample, NumStepSGD);
        const T result = dopant.template solve<RandomSource>(Target, params, dqmc);
        expect<RandomSource>(T(From) <= result && result <= T(To));

        expect(std::ranges::contains(dopant.getChemMus(), result));
        for (size_t i = 0; i < dopant.getNumSamples(); ++i) {
            if (result == dopant.getChemMus()[i]) {
                auto mean = dopant.getDensity()[i];
                auto devia = dopant.getNoises()[i];
                const T diff = abs(mean - Target);
                expect<RandomSource>(diff < T(0.01));
                expect<RandomSource>(diff < T(3) * devia);
                break;
            }
        }
    }
}

int main() {
    halfFilling();
    doping();
    return 0;
}
