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
using RandomSource = Random<>;
using T = float64;
using Tc = cfloat64;
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
    ActionMatrix<Tc> action(params, NumFreq, maxBoson);
    action.randAuxField<Random<>>();

    const MatrixND<Tc> result = action;
    size_t order = action.getOrder();
    for (size_t r = 0; r < order; ++r)
        for (size_t c = 0; c < order; ++c)
            expect(action.calc(r, c) == result[r, c]);
    return 0;
}