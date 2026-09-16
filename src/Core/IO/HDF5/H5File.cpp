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

namespace {
    class FileProperty {
        H5ID id;
    public:
        FileProperty();
        ~FileProperty() = default;
        FileProperty(const FileProperty&) = delete;
        FileProperty(FileProperty&&) noexcept = default;
        /* Getters */
        [[nodiscard]] auto getHID() const noexcept { return id.getHID(); }
    };
}

FileProperty::FileProperty() : id(H5Pcreate(H5P_FILE_ACCESS)) {
    std::ignore = H5Pset_fclose_degree(id.getHID(), H5F_CLOSE_STRONG);
}

H5File::H5File(H5ID id_) noexcept : H5Loc(std::move(id_)) {
    assert(Base::isa<H5File>());
}

bool H5File::isReadOnly() const noexcept {
    unsigned intent = H5F_ACC_RDWR;
    H5Fget_intent(getHID(), &intent);
    return (intent & H5F_ACC_RDWR) == 0;
}

H5File H5File::open(const char* name, unsigned int openflag) {
    FileProperty fapl;
    if (std::filesystem::exists(name)) {
        if (openflag & Trunc)
            return H5File(H5ID(H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl.getHID())));
        unsigned int access = (openflag & ReadWrite) ? H5F_ACC_RDWR : H5F_ACC_RDONLY;
        return H5File(H5ID(H5Fopen(name, access, fapl.getHID())));
    }
    if (!bool(openflag & ReadWrite))
        throw IOException("File not found");
    return H5File(H5ID(H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl.getHID())));
}
