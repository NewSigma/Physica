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
#include <filesystem>
#include "Physica/Core/IO/HDF5/H5File.h"
#include "Physica/Core/Exception/IOException.h"

using namespace Physica;

H5File::H5File(H5ID id_, unsigned int flag_) : H5Loc(std::move(id_)), openflag(flag_) {}

H5File H5File::open(const char* name, unsigned int openflag) {
    if (std::filesystem::exists(name)) {
        if (openflag & Trunc) {
            auto fid = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            return H5File(H5ID(fid), openflag);
        }
        unsigned int access = (openflag & ReadWrite) ? H5F_ACC_RDWR : H5F_ACC_RDONLY;
        auto fid = H5Fopen(name, access, H5P_DEFAULT);
        return H5File(H5ID(fid), openflag);
    }
    if (!bool(openflag & ReadWrite))
        throw IOException("File not found");
    auto fid = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    return H5File(H5ID(fid), openflag | Creat);
}
