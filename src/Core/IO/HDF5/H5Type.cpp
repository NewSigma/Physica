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
#include "Physica/Core/IO/HDF5/H5Type.h"
#include <cassert>

using namespace Physica;

H5Type::H5Type(hid_t hid) : H5ID(hid) {}

void H5Type::insert(const char* name, size_t offset, const H5Type& memberType) noexcept {
    assert(isCompound());
    H5Tinsert(getHID(), name, offset, memberType.getHID());
}

bool H5Type::isCompound() const noexcept {
    return H5Tget_class(getHID()) == H5T_COMPOUND;
}
