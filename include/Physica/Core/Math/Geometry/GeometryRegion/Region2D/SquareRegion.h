/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
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
