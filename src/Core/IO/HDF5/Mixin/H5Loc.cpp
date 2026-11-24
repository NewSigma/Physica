/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <utility>
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/IO/HDF5/H5Group.h"

using namespace Physica;

H5Loc::H5Loc(H5ID id_) noexcept : H5ID(std::move(id_)) {}

bool H5Loc::exists(const char* name) const {
    return H5Lexists(getHID(), name, H5P_DEFAULT) > 0;
}

H5Group H5Loc::openGroup(const char* name) {
    if (exists(name))
        return H5Group::open(*this, name);
    if (isReadOnly())
        throw IOException("[Error]: Group not found");
    return H5Group::create(*this, name);
}

H5Group H5Loc::openGroup(const char* name) const {
    if (!exists(name))
        throw IOException("[Error]: Group not found");
    return H5Group::open(*this, name);
}
