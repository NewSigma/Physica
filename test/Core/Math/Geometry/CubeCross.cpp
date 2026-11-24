/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Geometry/CubeCross.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;

int main() {
    {
        const ScalarType rep3 = reciprocal(ScalarType(3));
        const auto cross = CubeCross<ScalarType>({rep3, rep3, rep3, 1});
        expect(scalarNear(cross.getArea(), ScalarType(0), 1E-15));
    }
    {
        const auto cross = CubeCross<ScalarType>({2, 0, 0, 1});
        expect(scalarNear(cross.getArea(), ScalarType(4), 1E-15));
    }
    {
        const auto cross = CubeCross<ScalarType>({-0.5, 0, 0, 1});
        expect(scalarNear(cross.getArea(), ScalarType(0), 1E-15));
    }
    {
        const auto cross = CubeCross<ScalarType>({1, 1, 0, 1});
        expect(scalarNear(cross.getArea(), ScalarType(2 * M_SQRT2), 1E-15));
    }
    {
        const auto cross = CubeCross<ScalarType>({1, 1, 1, 1});
        expect(scalarNear(cross.getArea(), sqrt(ScalarType(3)) * ScalarType(2), 1E-15));
    }
    {
        const ScalarType rep3 = reciprocal(ScalarType(3));
        const ScalarType sqrt3 = sqrt(ScalarType(3));
        const auto x = VectorND<ScalarType>::linspace(-1, 1, 200);
        VectorND<ScalarType> y(x.getLength());
        for (size_t i = 0; i < y.getLength(); ++i) {
            Vector4D<ScalarType> plane{rep3, rep3, rep3, x[i]};
            plane *= ScalarType(10);
            const auto cross = CubeCross<ScalarType>(plane);
            y[i] = cross.getArea() * sqrt3;
        }
        const ScalarType volume = y.sum() * (x[1] - x[0]);
        expect(scalarNear(volume, ScalarType(8), 1E-6));
    }
    return 0;
}
