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
#include "Physica/Core/IO/HDF5/H5Attribute.h"
#include <cassert>
#include <utility>
#include "H5Apublic.h"
#include "H5Spublic.h"

using namespace Physica;

H5Attribute::H5Attribute(H5ID id) noexcept : H5ID(std::move(id)) {
    assert(isa<H5Attribute>());
}

H5Attribute& H5Attribute::operator=(H5Attribute obj) noexcept {
    swap(obj);
    return *this;
}

void H5Attribute::read(const H5Type& dtype, void* buf) const {
    H5Aread(getHID(), dtype.getHID(), buf);
}

void H5Attribute::write(const H5Type& dtype, const void* buf) const {
    checkWrite();
    H5Awrite(getHID(), dtype.getHID(), buf);
}

size_t H5Attribute::getSize() const noexcept {
    const auto space = H5ID(H5Aget_space(getHID()));
    const auto dim = H5Sget_simple_extent_ndims(space.getHID());
    assert(dim >= 0 && "[Error]: Failed to query the attribute space");
    assert(dim <= 1 && "[Error]: High dim attribute is not supported");
    hsize_t result = 1;
    if (dim == 1)
        H5Sget_simple_extent_dims(space.getHID(), &result, nullptr);
    return result;
}

void H5Attribute::swap(H5Attribute& obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    H5ID::swap(obj);
}
