/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
