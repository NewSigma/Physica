/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/BasicAlgebra/Polynomial.h"
#include "Test.h"

using namespace Physica;

namespace {
    template<Vector V>
    void testPolyRoot(const V& coeffs, double precision) {
        using ScalarType = V::ScalarType;
        const Polynomial<ScalarType, coeffs.getSizeAtCompile()> poly(coeffs);
        auto roots = polyRoot(poly);
        for (const auto& root : roots) {
            auto result = poly(root);
            expect(scalarNear(result, decltype(result)(0), precision));
        }
    }
}

int main() {
    using ScalarType = float64;
    VectorND<ScalarType> coeffs{1, 2, 3, 4, 5, 6};
    testPolyRoot(coeffs, 1E-11);
    return 0;
}
