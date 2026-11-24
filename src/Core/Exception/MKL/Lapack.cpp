/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <format>
#include <system_error>
#include "../SystemError.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

using namespace Physica;

namespace {
    class Impl final : public std::error_category {
    public:
        Impl() = default;
        Impl(const Impl&) = delete;
        Impl(Impl&&) noexcept = delete;
        ~Impl() = default;
        /* Operators */
        Impl& operator=(const Impl&) = delete;
        Impl& operator=(Impl&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] const char* name() const noexcept final { return "Intel MKL Lapack"; }
        [[nodiscard]] std::string message(int code) const final { return std::format("error code: {}", code); }
    };
}

void Physica::check_lapack(int err) {
    constinit static const Impl instance{};
    if (err != 0) [[unlikely]]
        throw SystemError(err, instance);
}
