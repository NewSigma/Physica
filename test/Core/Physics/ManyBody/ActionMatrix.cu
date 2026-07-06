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
#include "Physica/Core/Physics/ManyBody/DQMCImpl/ActionMatrix.cuh"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;
using T = float32;
using Tc = cfloat32;
constexpr int Dim = 2;
constexpr T HoppingT = 1;
constexpr T RepelU = 8;
constexpr T Beta = 8;
constexpr T ChemMu = -2;
constexpr int NumSplit = 1;

namespace {
    void hostDeviceMatch() {
        constexpr Array<int, 3> SmallInts{1, 2, 3};
        const size_t numSite = SmallInts.template select<RandomSource>() + 1;
        const int numFreq = SmallInts.template select<RandomSource>();

        const SquareLattice<Dim> lattice({numSite, numSite}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, ChemMu, NumSplit);
        int maxBoson = std::uniform_int_distribution<>(1, numFreq)(RandomSource::getInstance());
        device_obj<ActionMatrix<Tc>> d_action(params, numFreq, maxBoson);
        d_action.random_normal<Random<>>();
        device_obj<MatrixND<Tc>> result = d_action;

        ActionMatrix<Tc> action(params, numFreq, maxBoson);
        d_action.getAuxField().toHostAsync(action.getAuxField());
        MatrixND<Tc> answer = action;

        expect(matrixNear(answer, result.toHost(), 1E-6));
    }
}

int main() {
    hostDeviceMatch();
    return 0;
}
