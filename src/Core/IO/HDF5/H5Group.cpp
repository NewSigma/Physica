/*
 * Copyright 2026 Weibo He.
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
#include <cassert>
#include "Physica/Core/IO/HDF5/H5Group.h"

using namespace Physica;

H5Group::H5Group(H5ID id) : H5Loc(std::move(id)) {}

H5Group H5Group::create(const H5ID& loc, const char* name) {
    assert(!loc.isReadOnly());
    return H5Group(H5ID(H5Gcreate2(loc.getHID(), name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)));
}

H5Group H5Group::open(const H5ID& loc, const char* name) {
    return H5Group(H5ID(H5Gopen2(loc.getHID(), name, H5P_DEFAULT)));
}
