/*
 * Copyright 2022 WeiBo He.
 *
 * This file is part of Physica.

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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Geometry/GeoBase.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Float, false>;
using GeoBase2D = GeoBase<2, ScalarType>;

int main() {
    constexpr float epsilon = std::numeric_limits<float>::epsilon();
    if (!scalarNear(GeoBase2D::distToSegment({-1, -1}, {0, 0}, {1, 0}), sqrt(ScalarType::Two()), epsilon))
        return 1;
    if (!scalarNear(GeoBase2D::distToSegment({2, 1}, {0, 0}, {1, 0}), sqrt(ScalarType::Two()), epsilon))
        return 1;
    if (!scalarNear(GeoBase2D::distToSegment({1, -1}, {0, 0}, {1, 0}), ScalarType::One(), epsilon))
        return 1;
    return 0;
}
