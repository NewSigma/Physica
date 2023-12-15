/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    class H5File : public H5::H5File, public H5Location {
        using Base = H5::H5File;
        using Location = Core::H5Location;
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
    public:
        H5File(const char* name,
               unsigned int openflag,
               const H5::FileCreatPropList& create_plist = H5::FileCreatPropList::DEFAULT,
               const H5::FileAccPropList& access_plist = H5::FileAccPropList::DEFAULT);
        H5File(const H5File& obj);
        H5File(H5File&&) noexcept = delete;
        virtual ~H5File() = default;
        /* Operators */
        H5File& operator=(H5File& obj);
        H5File& operator=(H5File&&) noexcept = delete;
        /* Operations */
        using Location::exists;
        using Location::createDataSet;
        using Location::openDataSet;
        H5DataSet<1> createDataSet(const char* filepath, const char* name);
    };
}
