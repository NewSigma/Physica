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
#include "Physica/Core/Physics/ManyBody/FreqDQMC.cuh"
#include "Physica/Core/Physics/MC/HamiltonMC.h"
#include "Test.h"

using namespace Physica;
using T = float32;
using Tc = cfloat32;
using RandomSource = Random<>;
constexpr T HoppingT = 1;
constexpr T RepelU = MathConst<T>::pi;
constexpr T Beta = MathConst<T>::e;
constexpr int Dim = 1;
constexpr T StepSize = 1E-3;
constexpr T Duration = 1; // Smaller than host as it too slow

namespace {
    void hostDeviceCross() {
        const size_t numSite = Array<int, 2>{2, 3}.select<RandomSource>();
        const int freqDensity = Array<int, 2>{0, 1}.select<RandomSource>();

        const SquareLattice<Dim> lattice({numSite}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, 1);
        auto dqmcH = FreqDQMC<Tc>(params, freqDensity);
        auto dqmcD = device_obj<FreqDQMC<Tc>>(params, freqDensity);
        dqmcH.step_random<RandomSource>();

        VectorND<T> pos(dqmcH.getAuxField().getSize() * 2);
        pos.read(dqmcH.getAuxField().flatten());
        {
            auto potH = dqmcH.potentialV(pos);
            auto potD = dqmcD.potentialV(pos);
            expect<RandomSource>(scalarNear(potH, potD, 1E-6));
        }
        {
            VectorND<T> forceH, forceD;
            forceH.resize(pos);
            forceD.resize(pos);
            dqmcH.forceAsync(pos, forceH);
            dqmcD.forceAsync(pos, forceD);
            expect<RandomSource>(vectorNear(forceH, forceD, 1E-3));
        }
    }

    void conserve() {
        const size_t numSite = Array<int, 2>{2, 3}.select<RandomSource>();
        const int freqDensity = Array<int, 2>{0, 1}.select<RandomSource>();

        const SquareLattice<Dim> lattice({numSite}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, 1);
        auto dqmc = device_obj<FreqDQMC<Tc>>(params, freqDensity);
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
        constexpr T HoppingT = 0;
        constexpr T RepelU = 0;
        const T beta = T::random_uniform<RandomSource>() * 8;
        const SquareLattice<2> lattice({2, 2}, 1);
        HubbardParams<T> params(HoppingT, RepelU, lattice, beta, RepelU * 0.5, 1);
        {
            auto dqmc = device_obj<FreqDQMC<Tc>>(params, 1);
            dqmc.step_random<RandomSource>();
            dqmc.step<RandomSource>();
            expect<RandomSource>(dqmc.calcBerry().isZero());
        }
        params.setChemMu(T::random_uniform<RandomSource>() * 8 - 4);
        {
            auto dqmc = device_obj<FreqDQMC<Tc>>(params, 1);
            dqmc.step_random<RandomSource>();
            dqmc.step<RandomSource>();
            expect<RandomSource>(scalarNear(dqmc.calcBerry(), T(0), 1E-5));
        }
    }
}

int main() {
    hostDeviceCross();
    conserve();
    berry();
    return 0;
}
