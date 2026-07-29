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
#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/Utils/Container/Array.h"

using namespace Physica;

H5File::H5File(H5ID id_, unsigned int flag_) : H5Loc(std::move(id_)), openflag(flag_) {}

H5DataSet<1> H5File::createDataSet(const std::string& filepath, const std::string& name) {
    std::ifstream fin(filepath);
    if (!fin)
        throw IOException("[Error]: File not found");
    fin.seekg(0, std::ios::end);
    const auto size = fin.tellg();
    fin.seekg(0, std::ios::beg);

    auto buffer = Array<char>(size);
    fin.read(buffer.data(), size);
    const auto space = H5DataSpace<1>(size);
    auto dataset = createDataSet<1>(name, H5Type::get<char>(), space);
    dataset.write(buffer.data(), H5Type::get<char>());
    return dataset;
}

H5Group H5File::openGroup(const std::string& name) {
    if (exists(name))
        return H5Group::open(*this, name);
    if (isReadOnly())
        throw IOException("[Error]: Group not found");
    return H5Group::create(*this, name);
}

H5File H5File::open(const std::string& name, unsigned int openflag) {
    if (std::filesystem::exists(name)) {
        if (openflag & Trunc) {
            auto fid = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            return H5File(H5ID(fid), openflag);
        }
        unsigned int access = (openflag & ReadWrite) ? H5F_ACC_RDWR : H5F_ACC_RDONLY;
        auto fid = H5Fopen(name.c_str(), access, H5P_DEFAULT);
        return H5File(H5ID(fid), openflag);
    }
    if (!bool(openflag & ReadWrite))
        throw IOException("File not found");
    auto fid = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    return H5File(H5ID(fid), openflag | Creat);
}
