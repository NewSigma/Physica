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
#include "Physica/Core/Physics/MC/HamiltonMC.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;
using RandomSource = Random<>;
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

        expect<RandomSource>(scalarNear(prevE, curE, 1E-4)); // Energe conserves
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
    berry();
    return 0;
}
