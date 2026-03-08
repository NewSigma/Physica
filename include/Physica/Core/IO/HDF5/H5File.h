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
#pragma once

#include "HDF5.h"

namespace Physica {
    class PHYSICA_API H5File : public H5::H5File, public H5Loc {
        using This = H5File;
        using Base = H5::H5File;
        using Loc = H5Loc;
    public:
        enum OpenFlag : unsigned int {
            ReadOnly = 0x0000U,
            ReadWrite = 0x0001U,
            Trunc = 0x0002U,
            Excl = 0x0004U,
            Creat = 0x0010U,
            SingleWriteMultiRead_Write = 0x0020U,
            SingleWriteMultiRead_Read = 0x0040U
        };
    private:
        unsigned int openflag;
    public:
        H5File() = default;
        H5File(const This& obj) = delete;
        H5File(This&&) noexcept = delete;
        virtual ~H5File() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        using Loc::exists;
        using Loc::createDataSet;
        using Loc::openDataSet;
        using Loc::openGroup;
        [[nodiscard]] H5DataSet<1> createDataSet(const std::string& filepath, const std::string& name);
        [[nodiscard]] H5Group openGroup(const std::string& name);
        /* Getters */
        [[nodiscard]] unsigned int getOpenflag() const noexcept { return openflag; }
        [[nodiscard]] bool isReadOnly() const noexcept { return (openflag & ReadWrite) == 0; }
        /* Static members */
        [[nodiscard]] static H5File open(const std::string& name, unsigned int openflag = ReadWrite);
    private:
        H5File(const std::string& name,
               unsigned int openflag_ = OpenFlag::ReadOnly,
               const H5::FileCreatPropList& create_plist = H5::FileCreatPropList::DEFAULT,
               const H5::FileAccPropList& access_plist = H5::FileAccPropList::DEFAULT);
    };
}
