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
#include "Physica/Core/Physics/ManyBody/DQMCImpl/ActionMatrix.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<PCG64DXSM, 1000>;
using T = float64;
using Tc = cfloat64;
constexpr int Dim = 2;
constexpr T HoppingT = 1;
constexpr T RepelU = 8;
constexpr T Beta = 8;
constexpr T ChemMu = -2;
constexpr Array<int, 3> SmallInts{1, 2, 3};
constexpr int NumSplit = 1;

namespace {
    void assign_calc_match() {
        const size_t numSite = SmallInts.template select<RandomSource>() + 1;
        const int numFreq = SmallInts.template select<RandomSource>();
        const SquareLattice<Dim> lattice({numSite, numSite}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, ChemMu, NumSplit);
        int maxBoson = std::uniform_int_distribution<>(1, numFreq)(RandomSource::getInstance());
        ActionMatrix<Tc> action(params, numFreq, maxBoson);
        action.randAuxField<RandomSource>();

        const MatrixND<Tc> result = action;
        size_t order = action.getOrder();
        for (size_t r = 0; r < order; ++r)
            for (size_t c = 0; c < order; ++c)
                expect(action.calc(r, c) == result[r, c]);
    }

    void assign_gemv_match() {
        const size_t numSite = SmallInts.template select<RandomSource>() + 1;
        const int numFreq = SmallInts.template select<RandomSource>();
        const SquareLattice<Dim> lattice({numSite, numSite}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, ChemMu, NumSplit);
        int maxBoson = std::uniform_int_distribution<>(1, numFreq)(RandomSource::getInstance());
        ActionMatrix<Tc> action(params, numFreq, maxBoson);
        action.randAuxField<RandomSource>();

        size_t order = action.getOrder();
        const MatrixND<Tc> answer = action;
        auto result = VectorND<Tc>(order);
        for (size_t i = 0; i < order; ++i) {
            result = action * UnitVector<T>(i, order);
            expect(vectorNear(answer.col(i), result, 2UL));
        }
    }
}

int main() {
    assign_calc_match();
    assign_gemv_match();
    return 0;
}