/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
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