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
#include "Physica/Core/IO/HDF5/Mixin/H5ID.h"
#include <cassert>

using namespace Physica;

H5ID::H5ID(hid_t id) noexcept : id(id) {}

H5ID::H5ID(const H5ID& other) noexcept : id(other.id) {
    if (isValid())
        incRef();
}

H5ID::H5ID(H5ID&& other) noexcept : id(std::exchange(other.id, H5I_INVALID_HID)) {}

H5ID::~H5ID() noexcept {
    if (isValid())
        H5Idec_ref(id);
}

H5ID& H5ID::operator=(H5ID other) noexcept {
    swap(other);
    return *this;
}

bool H5ID::isValid() const noexcept {
    return H5Iis_valid(id) > 0;
}

bool H5ID::isFile() const noexcept {
    return H5Iget_type(id) == H5I_FILE;
}

bool H5ID::isGroup() const noexcept {
    return H5Iget_type(id) == H5I_GROUP;
}

bool H5ID::isDatatype() const noexcept {
    return H5Iget_type(id) == H5I_DATATYPE;
}

bool H5ID::isDataspace() const noexcept {
    return H5Iget_type(id) == H5I_DATASPACE;
}

bool H5ID::isDataset() const noexcept {
    return H5Iget_type(id) == H5I_DATASET;
}

bool H5ID::isAttribute() const noexcept {
    return H5Iget_type(id) == H5I_ATTR;
}

void H5ID::incRef() const noexcept {
    assert(isValid());
    [[maybe_unused]] auto err = H5Iinc_ref(id);
    assert(err >= 0);
}

void H5ID::swap(H5ID& other) noexcept {
    assert(this != &other && "[Error]: Self swap is likely a bug");
    std::swap(id, other.id);
}
