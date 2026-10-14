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
#include "Physica/Core/Physics/ManyBody/FreqDQMC.h"
#include "Physica/Core/Physics/ManyBody/ElasticDQMC.h"
#include "Physica/Core/Physics/ManyBody/GreenSampler/ScalarSampler.h"
#include "Physica/Core/Physics/MC/HamiltonMC.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;
using RandomSource = Random<PCG64DXSM, 566797295274468060>;
constexpr int Dim = 1;
constexpr T StepSize = 1E-3;
constexpr T Duration = 10;

namespace {
    void conserve() {
        constexpr T HoppingT = 1;
        constexpr T RepelU = MathConst<T>::pi;
        constexpr T Beta = MathConst<T>::e;
        const size_t numSite = Array<int, 2>{2, 3}.select<RandomSource>();
        const int freqDensity = Array<int, 2>{0, 1}.select<RandomSource>();
        const SquareLattice<Dim> lattice({numSite}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, 1);
        auto dqmc = FreqDQMC<Tc>(params, freqDensity);
        dqmc.step_for<RandomSource>(1);

        auto& engine = dqmc.getHMC().getRoot();
        engine.setTimeStep(StepSize);

        using Kinetic = OpenModel<T, 1, 1>;
        auto kinetic = Kinetic(1, 1);
        engine.template initMomentum<Kinetic, RandomSource>();

        const T prevE = engine.calcClassicalInternalEnergy(dqmc);
        engine.nve_step_for(Duration, kinetic, dqmc);
        const T curE = engine.calcClassicalInternalEnergy(dqmc);

        expect<RandomSource>(scalarNear(prevE, curE, 2E-4)); // Energy conserves
    }

    void elastic() {
        // ElasticDQMC is exact for the U = 0
        constexpr T HoppingT = 1;
        constexpr T RepelU = 0;
        constexpr T Beta = 4;
        constexpr int NumSite = 3;
        constexpr int FreqDensity = 1;
        const T chemMu = T::random_uniform<RandomSource>() * 2 - 1;

        const SquareLattice<Dim> lattice({NumSite}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, chemMu, 1);
        auto dqmc = FreqDQMC<Tc>(params, FreqDensity);
        dqmc.step_random<RandomSource>();
        dqmc.step<RandomSource>();

        ElasticDQMC<T> exact(params, FreqDensity);
        exact.step_random<RandomSource>();
        for (int spin : {0, 1})
            expect<RandomSource>(matrixNear(dqmc.getGreens()[spin], exact.getGreens()[spin], 1E-10));
    }

    void green() {
        // Regression test for green's function; verify we fo not mix non-diagonal terms into the equal-time green's function
        constexpr T HoppingT = 1;
        constexpr T RepelU = 2;
        constexpr T Beta = 4;
        constexpr T ChemMu = RepelU * 0.5;
        constexpr int NumSiteX = 2;
        constexpr int NumSiteY = 2;
        constexpr int FreqDensity = 1;
        constexpr int NumWarmup = 256;
        constexpr int NumSample = NumWarmup * 16;
        constexpr double Prec = 2E-2;

        const SquareLattice<2> lattice({NumSiteX, NumSiteY}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, ChemMu, 1);
        auto calcDoubleOccupy = [&](int maxBoson) {
            auto dqmc = FreqDQMC<Tc>(params, FreqDensity, maxBoson);
            dqmc.step_random<RandomSource>();
            dqmc.step_for<RandomSource>(NumWarmup);

            ScalarSampler<T> density(params, NumSample), doubleOccupy(params, NumSample);
            for (int i = 0; i < NumSample; ++i) {
                dqmc.step<RandomSource>();
                density.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<T>::Density);
                doubleOccupy.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<T>::DoubleOccupy);
            }
            expect<RandomSource>(scalarNear(density.calcMean(), T(1), 1E-10));
            return doubleOccupy.calcMean();
        };
        const T elastic = calcDoubleOccupy(1);
        const T full = calcDoubleOccupy(2 * ElasticDQMC<T>::calcFreqCutoff(Beta, FreqDensity));
        expect<RandomSource>(scalarNear(elastic, full, Prec));
    }

    void berry() {
        const SquareLattice<Dim> lattice({2}, 1);
        {
            constexpr T HoppingT = 0;
            constexpr T RepelU = 0;
            const T beta = T::random_uniform<RandomSource>() * 8;
            const T chemMu = 0;
            HubbardParams<T> params(HoppingT, RepelU, lattice, beta, chemMu, 1);
            auto dqmc = FreqDQMC<Tc>(params, 1);
            dqmc.step_random<RandomSource>();
            dqmc.step<RandomSource>();
            expect<RandomSource>(dqmc.calcBerry().isZero()); // Berry is trivially zero if t = U = \mu = 0

            params.setChemMu(T::random_uniform<RandomSource>() * 8 - 4);
            dqmc.step_random<RandomSource>();
            dqmc.step<RandomSource>();
            expect<RandomSource>(scalarNear(dqmc.calcBerry(), T(0), 1E-14));
        }
        {
            constexpr T HoppingT = 1;
            constexpr T RepelU = 1;
            const T beta = 0;
            const T chemMu = T::random_uniform<RandomSource>() * 8 - 4;
            HubbardParams<T> params(HoppingT, RepelU, lattice, beta, RepelU * 0.5 + chemMu, 1);
            auto dqmc = FreqDQMC<Tc>(params, 1);
            dqmc.step_random<RandomSource>();
            dqmc.step<RandomSource>();
            expect<RandomSource>(dqmc.calcBerry().isZero()); // Berry is trivially zero if beta is zero
        }
    }
}

int main() {
    conserve();
    elastic();
    green();
    berry();
    return 0;
}
