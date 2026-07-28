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

#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Core/Exception/IOException.h"

namespace Physica {
    class H5Group;

    class PHYSICA_API H5Loc {
    public:
        ~H5Loc() = default;
        /* Operations */
        [[nodiscard]] inline bool exists(const std::string& name, const H5::LinkAccPropList& lapl = H5::LinkAccPropList::DEFAULT) const;

        template<size_t Dim>
        H5DataSet<Dim> createDataSet(
                const std::string&name,
                const H5::DataType& data_type,
                const H5::DataSpace& data_space,
                const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT,
                const H5::DSetAccPropList& dapl = H5::DSetAccPropList::DEFAULT,
                const H5::LinkCreatPropList& lcpl = H5::LinkCreatPropList::DEFAULT) const;

        template<size_t Dim>
        [[nodiscard]] H5DataSet<Dim> openDataSet(const std::string& name);
        template<size_t Dim>
        [[nodiscard]] const H5DataSet<Dim> openDataSet(const std::string& name) const;

        [[nodiscard]] const H5Group openGroup(const std::string& name) const;
        [[nodiscard]] H5Group openGroup(const std::string& name);
    protected:
        H5Loc() = default;
        H5Loc(const H5Loc&) = default;
        H5Loc(H5Loc&&) noexcept = default;
        /* Operators */
        H5Loc& operator=(const H5Loc&) = default;
        H5Loc& operator=(H5Loc&&) noexcept = default;
    private:
        [[nodiscard]] auto& getDerived() { return *reinterpret_cast<H5::H5Location*>(this); }
        [[nodiscard]] const auto& getDerived() const { return *reinterpret_cast<const H5::H5Location*>(this); }
    };

    inline bool H5Loc::exists(const std::string& name, const H5::LinkAccPropList& lapl) const {
        return getDerived().exists(name, lapl);
    }

    template<size_t Dim>
    H5DataSet<Dim> H5Loc::createDataSet(
            const std::string&name,
            const H5::DataType& data_type,
            const H5::DataSpace& data_space,
            const H5::DSetCreatPropList& create_plist,
            const H5::DSetAccPropList& dapl,
            const H5::LinkCreatPropList& lcpl) const {
        return getDerived().createDataSet(name, data_type, data_space, create_plist, dapl, lcpl);
    }

    template<size_t Dim>
    H5DataSet<Dim> H5Loc::openDataSet(const std::string& name) {
        if (!exists(name))
            throw IOException("[Error]: Dataset not found");
        return getDerived().openDataSet(name);
    }

    template<size_t Dim>
    const H5DataSet<Dim> H5Loc::openDataSet(const std::string& name) const {
        if (!exists(name))
            throw IOException("[Error]: Dataset not found");
        return getDerived().openDataSet(name);
    }
}
