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
            expect<RandomSource>(scalarNear(mean, y_elem, 1E-12));
            expect<RandomSource>(devia.isZero());
        }
    }
}

int main() {
    samples();
    return 0;
}
