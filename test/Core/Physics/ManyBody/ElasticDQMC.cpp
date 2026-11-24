/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
