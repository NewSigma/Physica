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
#include <cmath>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/ForwardDenseQR.h"
#include "Physica/Core/Math/Transform/DiffFFT.h"
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
    void halfFillTest() {
        // Fuzzing test that half filling is free of sign problem
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

    void freeFermion() {
        // Test that trotter decomposition does not affect results of free system
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

    void complex() {
        // Test that particle number is real under complex DQMC
        constexpr double HoppingT = 1;
        constexpr double RepelU = 8;
        constexpr double Beta = 4;
        constexpr int NumSiteX = 4;
        constexpr int NumSiteY = 4;
        constexpr int NumSplit = Beta * 8;
        constexpr int NumSample = 64;
        const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
        const HubbardParams<Tc> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, NumSplit);
        auto dqmc = DQMC<Tc>(params);
        dqmc.step_random<RandomSource>();
        dqmc.step_for<RandomSource>(NumSample);
        for (const auto& green : dqmc.getGreens())
            for (auto elem : green.diag().view())
                expect<RandomSource>(scalarNear(elem.imag(), T(0), 1E-8));
    }

    void stability() {
        // Test that the Green's function's values have a manageable dynamic range
        constexpr static T StabilityMagnitudeLimit = T(1E100);
        constexpr int NumSite = 8;
        constexpr T RepelU = 10;
        constexpr T Beta = 60;

        int numSplit = std::max(int((Beta * 8 + T(0.5)).toMachine()), 2);
        if (numSplit % 2 != 0)
            numSplit += 1;

        const SquareLattice<Dim> lattice({size_t(NumSite), size_t(NumSite)}, 1);
        const HubbardParams<T> params(1, RepelU, lattice, Beta, RepelU * T(0.5), numSplit);
        auto dqmc = DQMC<T>(params);
        dqmc.step_random<RandomSource>();
        dqmc.step_for<RandomSource>(numSplit);
        expect<RandomSource>(std::ranges::all_of(dqmc.getGreens(), [](const auto& green) {
            return green.isFinite() && (abs_elem(green).max() < StabilityMagnitudeLimit);
        }));
    }

    void forward() {
        // Test differentiable samplers against the exact atomic limit,
        // where the Trotter decomposition is exact and the single-site solution is analytic
        using dfloat = Diff<T, DiffMode::Forward>;
        constexpr T RepelU = 2;
        constexpr T Beta = 2;
        constexpr T ChemMu = 1.5;
        constexpr int NumSite = 2;
        constexpr int NumSplit = 2;
        constexpr int NumWarmup = 256;
        constexpr int NumSample = NumWarmup * 16;
        constexpr double ValuePrec = 5E-2;
        constexpr double GradPrec = 5E-2;

        const HubbardParams<dfloat> params(MatrixND<T>::zeros(NumSite), dfloat(RepelU, 1), Beta, ChemMu, NumSplit);
        auto dqmc = DQMC<dfloat>(params);
        dqmc.step_random<RandomSource>();
        dqmc.step_for<RandomSource>(NumWarmup);

        ScalarSampler<dfloat> densitySampler(params, NumSample);
        ScalarSampler<dfloat> doubleSampler(params, NumSample);
        for (int i = 0; i < NumSample; ++i) {
            dqmc.step<RandomSource>();
            densitySampler.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<dfloat>::Density);
            doubleSampler.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<dfloat>::DoubleOccupy);
        }

        const T e1 = exp(Beta * ChemMu);
        const T e2 = exp(Beta * (T(2) * ChemMu - RepelU));
        const T z = T(1) + T(2) * e1 + e2;
        const T density = T(2) * (e1 + e2) / z;
        const T ddensity = -T(2) * Beta * e2 * (T(1) + e1) / square(z);
        const T doubleOccupy = e2 / z;
        const T ddoubleOccupy = -Beta * e2 * (T(1) + T(2) * e1) / square(z);

        const dfloat sign = densitySampler.calcRSign();
        expect(sign.value() == T(1));
        expect<RandomSource>(scalarNear(sign.grad(), T(0), 1E-8));

        const dfloat rawRho = densitySampler.calcRawMean();
        const dfloat meanRho = densitySampler.calcMean();
        expect<RandomSource>(scalarNear(rawRho.value(), T(density), ValuePrec));
        expect<RandomSource>(scalarNear(rawRho.grad(), T(ddensity), GradPrec));
        expect<RandomSource>(scalarNear(meanRho.value(), T(density), ValuePrec));
        expect<RandomSource>(scalarNear(meanRho.grad(), T(ddensity), GradPrec));

        const dfloat rawOcc = doubleSampler.calcRawMean();
        const dfloat meanOcc = doubleSampler.calcMean();
        expect<RandomSource>(scalarNear(rawOcc.value(), T(doubleOccupy), ValuePrec));
        expect<RandomSource>(scalarNear(rawOcc.grad(), T(ddoubleOccupy), GradPrec));
        expect<RandomSource>(scalarNear(meanOcc.value(), T(doubleOccupy), ValuePrec));
        expect<RandomSource>(scalarNear(meanOcc.grad(), T(ddoubleOccupy), GradPrec));
    }
}

int main() {
    halfFillTest();
    freeFermion();
    complex();
    stability();
    forward();
    return 0;
}
