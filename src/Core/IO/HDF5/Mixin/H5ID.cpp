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
#include <utility>
#include <H5Fpublic.h>
#include "Physica/Core/IO/HDF5/H5File.h"
#include "Physica/Core/Utils/NoImpl.h"

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

void H5ID::swap(H5ID& other) noexcept {
    assert(this != &other && "[Error]: Self swap is likely a bug");
    std::swap(id, other.id);
}

bool H5ID::isValid() const noexcept {
    return H5Iis_valid(id) > 0;
}

bool H5ID::isReadOnly() const noexcept {
    assert(isValid());
    return H5ID(H5Iget_file_id(id)).cast<H5File>().isReadOnly();
}

void H5ID::incRef() const noexcept {
    assert(isValid());
    [[maybe_unused]] auto err = H5Iinc_ref(id);
    assert(err >= 0);
}

auto H5ID::itype() const noexcept -> IdentifierType {
    switch (H5Iget_type(id)) {
    case H5I_FILE:
        return IdentifierType::File;
    case H5I_GROUP:
        return IdentifierType::Group;
    case H5I_DATATYPE:
        return IdentifierType::Datatype;
    case H5I_DATASPACE:
        return IdentifierType::Dataspace;
    case H5I_DATASET:
        return IdentifierType::Dataset;
    case H5I_ATTR:
        return IdentifierType::Attribute;
    default:
        return IdentifierType::Invalid;
    }
}
