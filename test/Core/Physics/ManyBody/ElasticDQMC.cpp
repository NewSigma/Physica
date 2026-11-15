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
#include "Physica/Core/Physics/ManyBody/ElasticDQMC.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<PCG64DXSM, 566797295274468060>;

namespace {
    void correction() {
        constexpr int NumTerm = 200000;
        for (const T beta : {T(0.5), T(4), MathConst<T>::pi}) {
            for (const T shift : {T(0.5), T(-0.5), T(2)}) {
                for (const int numFreq : {16, 64, 256}) {
                    const T c = beta * shift;
                    const T twoPi = T(2) * MathConst<T>::pi;
                    T omega = T(2 * numFreq + 1) * MathConst<T>::pi;
                    T exact = 0;
                    for (int i = 0; i < NumTerm; ++i) {
                        exact += reciprocal(square(c) + square(omega)) * T(-2) * c;
                        omega += twoPi;
                    }
                    const T approx = ElasticDQMC<T>::calcLocalCorrection(beta, T(0), shift, numFreq);
                    expect<RandomSource>(scalarNear(approx, exact, 2.0 / numFreq)); // The simple correction is O(1 / numFreq^2)
                }
            }
        }
    }
}

int main() {
    correction();
    return 0;
}
