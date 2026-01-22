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
constexpr int NumSiteX = 2;
constexpr int NumSiteY = 2;
constexpr int NumFreq = 2;
constexpr int NumSplit = 1;

int main() {
    const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
    const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, ChemMu, NumSplit);
    int maxBoson = std::uniform_int_distribution<>(1, NumFreq)(RandomSource::getInstance());
    device_obj<ActionMatrix<Tc>> d_action(params, NumFreq, maxBoson);
    d_action.randAuxField<Random<>>();
    device_obj<MatrixND<Tc>> result = d_action;

    ActionMatrix<Tc> action(params, NumFreq, maxBoson);
    d_action.getAuxField().toHostAsync(action.getAuxField());
    MatrixND<Tc> answer = action;

    expect(matrixNear(answer, result.toHost(), 1E-6));
    return 0;
}