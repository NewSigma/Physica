/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstdint>
#include <source_location>
#include <utility>
#include "Physica/Macro.h"
#include "Physica/Core/Math/Random/Random.h"

namespace Physica {
    PHYSICA_API void expect(bool pass, std::source_location loc = std::source_location::current(), uint64_t seed = 0) noexcept;

    template<RNG R>
    void expect(bool pass, std::source_location loc = std::source_location::current()) noexcept {
        expect(pass, std::move(loc), R::getInstance().getSeed());
    }

    consteval void syntax_only([[maybe_unused]] auto expr) noexcept {}
}
