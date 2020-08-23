/*
 * Copyright 2020 WeiBo He.
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
#include <utility>
#include "Physica/Core/Math/Geometry/GeometryRegion/Region2D/SquareRegion.h"

namespace Physica::Core {
    SquareRegion::SquareRegion(Point2D p, float length) : GeometryRegion<2>(Square), p(std::move(p)), length(length) {}

    SquareRegion::SquareRegion(const SquareRegion& region)
            : GeometryRegion<2>(Square), p(region.p), length(region.length) {}

    SquareRegion::SquareRegion(SquareRegion&& region) noexcept
            : GeometryRegion<2>(Square), p(std::move(region.p)), length(region.length) {}

    SquareRegion& SquareRegion::operator=(const SquareRegion& region) {
        if(this != &region) {
            p = region.p;
            length = region.length;
        }
        return *this;
    }

    SquareRegion& SquareRegion::operator=(SquareRegion&& region) noexcept {
        p = std::move(region.p);
        length = region.length;
        return *this;
    }

    GeometryRegion<2>* SquareRegion::release() {
        return new SquareRegion(std::move(*this));
    }
}