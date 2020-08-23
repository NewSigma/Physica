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
#ifndef PHYSICA_SQUAREREGION_H
#define PHYSICA_SQUAREREGION_H

#include "Physica/Core/Math/Geometry/Point.h"
#include "Physica/Core/Math/Geometry/GeometryRegion/GeometryRegion.h"

namespace Physica::Core {
    /*!
     * Imagine a rectangular coordinate system whose x-axis points to right and y-axis points to up.
     * Then @param p1 is the left bottom point of the square and @param length is the length of sides.
     */
    class SquareRegion : public GeometryRegion<2> {
        Point2D p;
        float length;
    public:
        SquareRegion(Point2D p, float length);
        SquareRegion(const SquareRegion& region);
        SquareRegion(SquareRegion&& region) noexcept;
        ~SquareRegion() override = default;
        /* Operators */
        SquareRegion& operator=(const SquareRegion& region);
        SquareRegion& operator=(SquareRegion&& region) noexcept;
    protected:
        GeometryRegion<2>* release() override;
    };
}

#endif
