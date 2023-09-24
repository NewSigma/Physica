/*
 * Copyright 2023 WeiBo He.
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

using namespace Physica::Core;
using ScalarType = Scalar<Double>;

int main() {
    {
        const ScalarType rep3 = reciprocal(ScalarType(3));
        const auto cross = CubeCross<ScalarType>({rep3, rep3, rep3});
        if (!scalarNear(cross.getArea(), ScalarType(0), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({2, 0, 0});
        if (!scalarNear(cross.getArea(), ScalarType(4), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({-0.5, 0, 0});
        if (!scalarNear(cross.getArea(), ScalarType(0), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({1, 1, 0});
        if (!scalarNear(cross.getArea(), ScalarType(2 * M_SQRT2), 1E-15))
            return 1;
    }
    {
        const auto cross = CubeCross<ScalarType>({1, 1, 1});
        if (!scalarNear(cross.getArea(), sqrt(ScalarType(3)) * ScalarType(2), 1E-15))
            return 1;
    }
    return 0;
}
