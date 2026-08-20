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

#include "H5Group.h"

namespace Physica {
    class PHYSICA_API H5File : public H5Loc {
        using This = H5File;
        using Base = H5Loc;
    public:
        enum OpenFlag : int8_t {
            ReadOnly = 0x0000U,
            ReadWrite = 0x0001U,
            Trunc = 0x0002U,
            Excl = 0x0004U,
            Creat = 0x0010U,
            SingleWriteMultiRead_Write = 0x0020U,
            SingleWriteMultiRead_Read = 0x0040U
        };
    public:
        H5File(const This& obj) = delete;
        H5File(This&&) noexcept = default;
        ~H5File() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = default;
        /* Getters */
        [[nodiscard]] bool isReadOnly() const noexcept;
        /* Static members */
        [[nodiscard]] static H5File open(const char* name, unsigned int openflag = ReadWrite);
        [[nodiscard]] constexpr static IdentifierType itype() noexcept { return IdentifierType::File; }
    private:
        H5File(H5ID id_) noexcept;

        friend class H5ID;
    };
}
