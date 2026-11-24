/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
        H5File(H5ID id_);

        friend class H5ID;
    };
}
