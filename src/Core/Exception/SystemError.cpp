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
