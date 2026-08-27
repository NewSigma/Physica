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

#include "Mixin/H5Loc.h"

namespace Physica {
    class PHYSICA_API H5Group : public H5Loc {
        using This = H5Group;
    public:
        H5Group() = default;
        explicit H5Group(H5ID id) noexcept;
        H5Group(const This&) = default;
        H5Group(This&&) noexcept = default;
        ~H5Group() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Static members */
        [[nodiscard]] static H5Group create(const H5Loc& loc, const char* name);
        [[nodiscard]] static H5Group open(const H5Loc& loc, const char* name);
        [[nodiscard]] constexpr static IdentifierType itype() noexcept { return IdentifierType::Group; }
    };
}
