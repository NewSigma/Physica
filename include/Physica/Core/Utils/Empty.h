/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

namespace Physica {
    /**
     * Typically used in std::conditional_t as the fallback type; everything regresses to no-op
     */
    class Empty {
        using This = Empty;
    public:
        constexpr Empty() = default;
        constexpr Empty(auto&&...) {}
        constexpr Empty(const This&) = default;
        constexpr Empty(This&&) noexcept = default;
        constexpr ~Empty() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        constexpr void swap(This&) noexcept {}
    };
}
