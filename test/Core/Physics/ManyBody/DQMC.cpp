/*
 * Copyright 2025 Weibo He.
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

using namespace Physica;
using T = float64;
using Tc = cfloat64;
using RandomSource = Random<>;
constexpr int Dim = 2;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr double Beta = 8;
constexpr int NumSiteX = 4;
constexpr int NumSiteY = 4;
constexpr int NumSplit = Beta * 8;
constexpr int NumSample = 1024;

namespace {
    /**
     * Test that half filling is free of sign problem
     */
    void halfFillTest() {
        const SquareLattice<Dim> lattice({NumSiteX, NumSiteY}, 1);
        const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, NumSplit);
        auto dqmc = DQMC<T>(params);
        dqmc.step_random<RandomSource>();
        for (int i = 0; i < NumSample; ++i) {
            dqmc.step_spin<RandomSource>();
            if (dqmc.getSign().isNegative())
                exit(EXIT_FAILURE);
        }
    }
}

int main() {
    halfFillTest();
    return 0;
}
