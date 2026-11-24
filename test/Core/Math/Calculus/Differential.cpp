/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Math/Calculus/SpecialFunctions.h"
#include "Test.h"

int main() {
    using namespace Physica;
    using ScalarType = float64;
    {
        auto func = [](const ScalarType& x) { return hermiteH(5, x); };
        auto result = Differential<ScalarType>::forward(func, ScalarType(3), ScalarType(1E-8));
        const ScalarType answer = ScalarType(8760);
        expect(scalarNear(result, answer, 1E-8));

        result = Differential<ScalarType>::backward(func, ScalarType(3), ScalarType(1E-8));
        expect(scalarNear(result, answer, 1E-7));

        result = Differential<ScalarType>::doublePoint(func, ScalarType(3), ScalarType(1E-5));
        expect(scalarNear(result, answer, 1E-10));

        result = Differential<ScalarType>::ridders(func, ScalarType(3), ScalarType(1.55));
        expect(scalarNear(result, answer, 1E-14));
    }
    return 0;
}
