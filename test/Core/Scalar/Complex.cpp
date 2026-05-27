/*
 * Copyright 2024-2026 Weibo He.
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
#include "Physica/Core/Scalar/Complex.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void simd() {
        auto x = SIMD<float64, 2>::random_uniform<Random<>>();
        auto y = SIMD<cfloat64, 2>::random_uniform<Random<>>();
        auto z = (x / y) * y;
        for (int i = 0; i < x.size(); ++i)
            expect<RandomSource>(scalarNear(x[i], z[i].real(), 2UL));
    }
}

int main() {
    static_assert(std::formattable<Complex<float64>, char>);
    simd();
    expect(ln1pexp(cfloat64(-1000, 0)).isZero()); // Test underflow
    return 0;
}
