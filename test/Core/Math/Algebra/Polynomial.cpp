/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/Polynomial.h"
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
