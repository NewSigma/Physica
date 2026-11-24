/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Geometry/GeometryRegion/GeometryRegion.h"

using namespace Physica;

GeometryRegion<2>::GeometryRegion(RegionType type) : type(type) {}
GeometryRegion<3>::GeometryRegion(RegionType type) : type(type) {}
