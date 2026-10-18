/*
 * Copyright 2026 Weibo He.
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

#include <cassert>
#include <ranges>
#include <utility>
#include <H5Ppublic.h>
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/IO/HDF5/H5Attribute.h"
#include "Physica/Core/IO/HDF5/H5DataSpace.h"

namespace Physica {
    class Attributable {
        using This = Attributable;
    public:
        ~Attributable() = default;
        /* Operations */
        [[nodiscard]] H5Attribute openAttribute(this auto&&, const char* name);
        [[nodiscard]] H5Attribute createAttribute(this const auto&, const char* name, const H5Type& dtype, const auto& space);

        void readAttr(this const auto&, const char* name, auto& value);
        void writeAttr(this auto&& self, const char* name, const auto& value);
        /* Getters */
        [[nodiscard]] bool attrExists(this const auto&, const char* name) noexcept;
    protected:
        Attributable() = default;
        Attributable(const This&) = default;
        Attributable(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    private:
        /* Static members */
        template<class T>
        constexpr static size_t calcNumElem(const T& value) noexcept;
    };

    H5Attribute Attributable::openAttribute(this auto&& self, const char* name) {
        H5ID id(H5Aopen(self.getHID(), name, H5P_DEFAULT));
        if (!id.isValid())
            throw IOException("[Error]: Attribute not found");
        return H5Attribute(std::move(id));
    }

    H5Attribute Attributable::createAttribute(this const auto& self, const char* name, const H5Type& dtype, const auto& space) {
        self.checkWrite();
        return H5Attribute(H5ID(H5Acreate2(self.getHID(), name, dtype.getHID(), space.getHID(), H5P_DEFAULT, H5P_DEFAULT)));
    }

    void Attributable::readAttr(this const auto& self, const char* name, auto& value) {
        self.openAttribute(name).read(value);
    }

    void Attributable::writeAttr(this auto&& self, const char* name, const auto& value) {
        using T = std::remove_cvref_t<decltype(value)>;
        const auto type = H5Type::get<T>();
        const auto space = H5DataSpace<1>(calcNumElem(value));
        H5Attribute attr;
        if (self.attrExists(name))
            attr = self.openAttribute(name);
        else
            attr = self.createAttribute(name, type, space);
        attr.write(value);
    }

    bool Attributable::attrExists(this const auto& self, const char* name) noexcept {
        return H5Aexists(self.getHID(), name) > 0;
    }

    template<class T>
    constexpr size_t Attributable::calcNumElem([[maybe_unused]] const T& value) noexcept {
        if constexpr (std::is_array_v<T>) {
            static_assert(std::rank_v<T> == 1, "[Error]: High dim array is not supported");
            return std::extent_v<T>;
        }
        else if constexpr (std::ranges::contiguous_range<T>)
            return std::ranges::size(value);
        else
            return 1;
    }
}
