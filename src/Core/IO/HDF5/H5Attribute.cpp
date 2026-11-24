/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
