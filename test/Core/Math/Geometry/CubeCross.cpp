/*
 * Copyright 2023 Weibo He.
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
#include <iostream>
#include "Physica/Core/Math/Geometry/CubeCross.h"

using namespace Physica;
using ScalarType = float64;

int main() {
    {
        const ScalarType rep3 = reciprocal(ScalarType(3));
        const auto cross = CubeCross<ScalarType>({rep3, rep3, rep3, 1});
        if (!scalarNear(cross.getArea(), ScalarType(0), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({2, 0, 0, 1});
        if (!scalarNear(cross.getArea(), ScalarType(4), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({-0.5, 0, 0, 1});
        if (!scalarNear(cross.getArea(), ScalarType(0), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({1, 1, 0, 1});
        if (!scalarNear(cross.getArea(), ScalarType(2 * M_SQRT2), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({1, 1, 1, 1});
        if (!scalarNear(cross.getArea(), sqrt(ScalarType(3)) * ScalarType(2), 1E-15))
            return 1;
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
        if (!scalarNear(volume, ScalarType(8), 1E-6))
            return 1;
    }
    return 0;
}
