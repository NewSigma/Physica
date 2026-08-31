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

#include <cstdint>
#include <source_location>
#include <utility>
#include "Physica/Macro.h"
#include "Physica/Core/Math/Random/Random.h"

namespace Physica {
    enum class Backend : char {
        Base,
        MKL
    };

    PHYSICA_API void expect(bool pass, std::source_location loc = std::source_location::current(), uint64_t seed = 0) noexcept;

    template<RNG R>
    void expect(bool pass, std::source_location loc = std::source_location::current()) noexcept {
        expect(pass, std::move(loc), R::getInstance().getSeed());
    }

    consteval void syntax_only([[maybe_unused]] auto expr) noexcept {}
}
