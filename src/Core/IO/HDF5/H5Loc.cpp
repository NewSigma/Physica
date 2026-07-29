/*
 * Copyright 2023-2026 Weibo He.
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
#include "Physica/Core/IO/HDF5/H5Group.h"

using namespace Physica;

H5Loc::H5Loc(H5ID id_) : H5ID(std::move(id_)) {}

bool H5Loc::exists(const std::string& name) const {
    return H5Lexists(getHID(), name.c_str(), H5P_DEFAULT) > 0;
}

H5Group H5Loc::openGroup(const std::string& name) {
    if (exists(name))
        return H5Group::open(*this, name);
    return H5Group::create(*this, name);
}

H5Group H5Loc::openGroup(const std::string& name) const {
    if (!exists(name))
        throw IOException("[Error]: Group not found");
    return H5Group::open(*this, name);
}
