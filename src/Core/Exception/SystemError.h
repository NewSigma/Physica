/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <string>
#include <system_error>

namespace Physica {
    class SystemError : public std::system_error {
        using This = SystemError;
        using Base = std::system_error;
        std::string buf;
    public:
        SystemError(int code, const std::error_category& category);
        SystemError(const This&) = default;
        SystemError(This&&) noexcept = default;
        ~SystemError() noexcept override = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Getters */
        [[nodiscard]] const char* what() const noexcept override;
    };
}
