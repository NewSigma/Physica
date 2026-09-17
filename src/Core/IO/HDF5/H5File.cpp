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
    auto hid = id.getHID();
    std::ignore = H5Pset_fclose_degree(hid, H5F_CLOSE_STRONG);
    if constexpr (HasMPI()) {
        // MPI might fork process and hold the lock; disable locking to avoid leaking.
        std::ignore = H5Pset_file_locking(hid, false, true);
    }
}

H5File::H5File(H5ID id_) : H5Loc(std::move(id_)) {
    if (!(Base::isValid() && Base::isa<H5File>()))
        throw IOException("[Error]: Failed to open HDF5 file; May be locked by another process?");
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
