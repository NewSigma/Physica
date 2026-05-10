/*
 * Copyright 2025-2026 Weibo He.
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

#include <cstdlib>
#include <print>
#include <source_location>
#include "Physica/Core/Math/Random/Random.h"

namespace Physica {
    [[noreturn]] inline void expect_fail(std::source_location loc, uint64_t seed) noexcept {
        std::println("Failed at file: {}:{}:{}", loc.file_name(), loc.line(), loc.column());
        std::println("          func: {}", loc.function_name());
        if (seed != 0)
            std::println("          seed: {}", seed);
        exit(EXIT_FAILURE);
    }

    inline void expect(bool pass, std::source_location loc = std::source_location::current()) noexcept {
        if (!pass) [[unlikely]]
            expect_fail(loc, 0);
    }

    template<RNG R>
    void expect(bool pass, std::source_location loc = std::source_location::current()) noexcept {
        if (!pass) [[unlikely]]
            expect_fail(loc, R::getInstance().getSeed());
    }

    consteval void syntax_only([[maybe_unused]] auto expr) noexcept {}
}
