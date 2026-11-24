/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/ML/GPModel.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;

namespace {
    void samples() {
        GPModel<T> gp(1);
        const auto x = VectorND<T>::linspace(-MathConst<T>::pi, MathConst<T>::pi, 8);
        const auto y = VectorND<T>(sin(x) + VectorND<T>::random_normal<RandomSource>(x.getLength()) * T(1E-3));
        // Test exact prediction
        gp.regression(x.transpose(), y, T(0));
        for (auto [x_elem, y_elem] : zip(x, y)) {
            auto [mean, devia] = gp.predict({x_elem});
            expect<RandomSource>(scalarNear(mean, y_elem, 1E-10));
            expect<RandomSource>(devia.isZero());
        }
    }
}

int main() {
    samples();
    return 0;
}
