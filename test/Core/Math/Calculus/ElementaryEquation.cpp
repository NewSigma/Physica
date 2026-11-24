/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Calculus/ElementaryEquation.h"
#include "Test.h"

int main() {
    using namespace Physica;
    using ScalarType = float64;
    auto func = [](const ScalarType& x) { return exp(x) - ScalarType(5); };
    {
        ScalarType result = bisectionMethod(func, ScalarType(0), ScalarType(3));
        const ScalarType answer = ScalarType(1.6094379124341003746);
        expect(scalarNear(result, answer, 1E-15));
    }
    {
        ScalarType result = secant(func, ScalarType(0), ScalarType(3), ScalarType(1E-9));
        const ScalarType answer = ScalarType(1.6094379124341003746);
        expect(scalarNear(result, answer, 1E-15));
    }
    return 0;
}
