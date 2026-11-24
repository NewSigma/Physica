/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Geometry/GeoBase2D.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float32;

int main() {
    constexpr float epsilon = std::numeric_limits<float>::epsilon();
    expect(scalarNear(GeoBase2D<ScalarType>::distToSegment({-1, -1}, {0, 0}, {1, 0}), sqrt(ScalarType(2)), epsilon));
    expect(scalarNear(GeoBase2D<ScalarType>::distToSegment({2, 1}, {0, 0}, {1, 0}), sqrt(ScalarType(2)), epsilon));
    expect(scalarNear(GeoBase2D<ScalarType>::distToSegment({1, -1}, {0, 0}, {1, 0}), ScalarType(1), epsilon));
    {
        using VectorType = Vector2D<ScalarType>;
        using Triangle = Array<VectorType, 3>;
        {
            Triangle poly{{0.57, 0}, {0.67, 0}, {0.57, 0.05}};
            bool flag1 = GeoBase2D<ScalarType>::pointOnPoly(VectorType{0, 0.3}, poly);
            bool flag2 = GeoBase2D<ScalarType>::pointOnPoly(VectorType{0.62, 0.01}, poly);
            expect(!flag1 && flag2);
        }
        {
            Triangle poly{{0.095, 0}, {0.095, 0.047}, {0, 0.047}};
            bool flag1 = GeoBase2D<ScalarType>::pointOnPoly(VectorType{0.36, 0}, poly);
            expect(!flag1);
        }
    }
    return 0;
}
