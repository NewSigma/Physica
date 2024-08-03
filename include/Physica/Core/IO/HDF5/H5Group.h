/*
 * Copyright 2023 Weibo He.
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
    class PHYSICA_API H5Group : public H5::Group, public H5Location {
        using Base = H5::Group;
        using Location = Core::H5Location;
    public:
        H5Group(H5::Group group) : Base(group) {}
        H5Group(const H5Group&) = default;
        H5Group(H5Group&&) noexcept = delete;
        ~H5Group() = default;
        /* Operators */
        H5Group& operator=(H5Group& obj) = default;
        H5Group& operator=(H5Group&&) noexcept = delete;
        /* Operations */
        using Location::createDataSet;
        using Location::openDataSet;

        inline H5Group createGroup(const char* name, size_t size_hint = 0) const;
        [[nodiscard]] inline H5Group openGroup(const char* name) const;
    private:
        using Base::createDataSet;
        using Base::openDataSet;
    };

    inline H5Group H5Group::createGroup(const char* name, size_t size_hint) const {
        return Base::createGroup(name, size_hint);
    }

    inline H5Group H5Group::openGroup(const char* name) const {
        return Base::openGroup(name);
    }
}
