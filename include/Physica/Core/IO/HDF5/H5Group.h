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

#include "H5Type.h"
#include "HDF5.h"

namespace Physica {
    class PHYSICA_API H5Group : public H5Loc {
        using This = H5Group;
        using Loc = H5Loc;
    public:
        H5Group() = default;
        explicit H5Group(H5ID id);
        H5Group(const This&) = default;
        H5Group(This&&) noexcept = default;
        ~H5Group() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        using Loc::exists;
        using Loc::createDataSet;
        using Loc::openDataSet;
        using Loc::openGroup;

        template<class T>
        void readAttr(const std::string& name, T& value) const;
        template<class T>
        void writeAttr(const std::string& name, T value);
        /* Static members */
        [[nodiscard]] static H5Group create(const H5ID& loc, const std::string& name);
        [[nodiscard]] static H5Group open(const H5ID& loc, const std::string& name);
    };

    template<class T>
    void H5Group::readAttr(const std::string& name, T& value) const {
        constexpr bool IsArray = std::is_array<T>::value;
        constexpr size_t NumElem = IsArray ? std::extent<T>::value : 1;
        static_assert(!IsArray || std::rank<T>::value == 1, "[Error]: High dim array is not supported");
        static_assert(NumElem > 0, "[Error]: Bad array size");

        const auto type = H5Type::get<T>();
        const auto space = H5DataSpace<1>(NumElem);
        H5Attribute attr;
        if (attrExists(name.c_str()))
            attr = openAttribute(name);
        else
            attr = createAttribute(name, type, space);
        attr.read(type, &value);
    }

    template<class T>
    void H5Group::writeAttr(const std::string& name, T value) {
        constexpr bool IsArray = std::is_array<T>::value;
        constexpr size_t NumElem = IsArray ? std::extent<T>::value : 1;
        static_assert(!IsArray || std::rank<T>::value == 1, "[Error]: High dim array is not supported");
        static_assert(NumElem > 0, "[Error]: Bad array size");

        const auto type = H5Type::get<T>();
        const auto space = H5DataSpace<1>(NumElem);
        H5Attribute attr;
        if (attrExists(name.c_str()))
            attr = openAttribute(name);
        else
            attr = createAttribute(name, type, space);
        attr.write(type, &value);
    }
}
