/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <cstddef>
#include <format>
#include <string>
#include <version>
#include "SystemError.h"

using namespace Physica;

namespace {
    std::string capture_stacktrace(std::size_t skip, std::size_t maxDepth) noexcept {
        // TODO: Wait for libc++ std::stacktrace
        std::ignore = skip;
        std::ignore = maxDepth;
        return {};
    }
}

SystemError::SystemError(int code, const std::error_category& category) : Base(code, category) {
    std::size_t skip = 1; // trace starts at the check function
    const auto trace = capture_stacktrace(skip, 32);
    if (!trace.empty())
        buf = std::format("{}\nStack trace:\n{}", Base::what(), trace);
}

const char* SystemError::what() const noexcept {
    return buf.empty() ? Base::what() : buf.c_str();
}
