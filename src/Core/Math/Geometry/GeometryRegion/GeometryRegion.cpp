/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <utility>
#include "Physica/Core/Math/Geometry/GeometryRegion/GeometryRegion.h"

namespace Physica::Core {
    /////////////////////////////////////Dimension 2/////////////////////////////////////
    GeometryRegion<2>::GeometryRegion(RegionType type) : type(type) {}
    /////////////////////////////////////Dimension 3/////////////////////////////////////
    GeometryRegion<3>::GeometryRegion(RegionType type) : type(type) {}
}