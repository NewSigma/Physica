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

#include <hdf5.h>
#include "H5Type.h"
#include "H5Attribute.h"
#include "Physica/Core/Exception/IOException.h"

namespace Physica {
    template<size_t Dim> class H5DataSpace;
    template<size_t Dim> class H5DataSet;
    class H5Group;

    class PHYSICA_API H5Loc : public H5ID {
        using This = H5Loc;
    public:
        H5Loc() = default;
        explicit H5Loc(H5ID id_);
        H5Loc(const This&) = default;
        H5Loc(This&&) noexcept = default;
        ~H5Loc() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<size_t Dim>
        [[nodiscard]] H5DataSet<Dim> createDataSet(const std::string& name, const H5Type& dtype, const H5DataSpace<Dim>& data_space) const;

        template<size_t Dim>
        [[nodiscard]] H5DataSet<Dim> openDataSet(const std::string& name);
        template<size_t Dim>
        [[nodiscard]] const H5DataSet<Dim> openDataSet(const std::string& name) const;

        [[nodiscard]] H5Group openGroup(const std::string& name);
        [[nodiscard]] H5Group openGroup(const std::string& name) const;

        [[nodiscard]] H5Attribute openAttribute(const std::string& name) const;
        template<class SpaceType>
        [[nodiscard]] H5Attribute createAttribute(const std::string& name, const H5Type& dtype, const SpaceType& space) const;
        /* Getters */
        [[nodiscard]] bool exists(const std::string& name) const;
        [[nodiscard]] bool attrExists(const char* name) const;
    };

    template<size_t Dim>
    H5DataSet<Dim> H5Loc::createDataSet(const std::string& name, const H5Type& dtype, const H5DataSpace<Dim>& data_space) const {
        return H5DataSet<Dim>(H5ID(H5Dcreate2(getHID(), name.c_str(), dtype.getHID(), data_space.getHID(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)));
    }

    template<size_t Dim>
    H5DataSet<Dim> H5Loc::openDataSet(const std::string& name) {
        if (!exists(name))
            throw IOException("[Error]: Dataset not found");
        return H5DataSet<Dim>(H5ID(H5Dopen2(getHID(), name.c_str(), H5P_DEFAULT)));
    }

    template<size_t Dim>
    const H5DataSet<Dim> H5Loc::openDataSet(const std::string& name) const {
        if (!exists(name))
            throw IOException("[Error]: Dataset not found");
        return H5DataSet<Dim>(H5ID(H5Dopen2(getHID(), name.c_str(), H5P_DEFAULT)));
    }

    template<class SpaceType>
    H5Attribute H5Loc::createAttribute(const std::string& name, const H5Type& dtype, const SpaceType& space) const {
        return H5Attribute(H5ID(H5Acreate2(getHID(), name.c_str(), dtype.getHID(), space.getHID(), H5P_DEFAULT, H5P_DEFAULT)));
    }
}
