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

#include "H5Type.h"

namespace Physica {
    class PHYSICA_API H5Attribute : public H5ID {
        using This = H5Attribute;
    public:
        H5Attribute() = default;
        explicit H5Attribute(H5ID id);
        H5Attribute(const This&) = default;
        H5Attribute(This&&) noexcept = default;
        ~H5Attribute() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void read(const H5Type& dtype, void* buf) const;
        void write(const H5Type& dtype, const void* buf) const;
        /* Getters */
        void swap(This& obj) noexcept;
    };
}
