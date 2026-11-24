/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <cassert>
#include <utility>
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/IO/HDF5/H5Group.h"

using namespace Physica;

H5Group::H5Group(H5ID id) noexcept : H5Loc(std::move(id)) {
    assert(isa<H5Group>());
}

H5Group H5Group::create(const H5Loc& loc, const char* name) {
    loc.checkWrite();
    return H5Group(H5ID(H5Gcreate2(loc.getHID(), name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)));
}

H5Group H5Group::open(const H5Loc& loc, const char* name) {
    H5ID id(H5Gopen2(loc.getHID(), name, H5P_DEFAULT));
    if (!id.isValid())
        throw IOException("[Error]: Group not found");
    return H5Group(std::move(id));
}
