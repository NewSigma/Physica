/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <array>
#include <cassert>
#include <cstdint>

namespace Physica {
    [[nodiscard]] constexpr uint64_t factorial(uint64_t x) noexcept {
        constexpr static std::array<uint64_t, 16> Table{
            1,
            1,
            2,
            6,
            24,
            120,
            720,
            5040,
            40320,
            362880,
            3628800,
            39916800,
            479001600,
            6227020800,
            87178291200,
            1307674368000,
        };

        assert(x <= 20 && "[Error]: factorial overflows");
        if (x < Table.size())
            return Table[x];

        uint64_t result = Table.back();
        for (uint64_t i = Table.size(); i <= x; ++i)
            result *= i;
        return result;
    }

    [[nodiscard]] constexpr uint64_t doubleFactorial(uint64_t x) noexcept {
        constexpr static std::array<uint64_t, 16> Table{
            1,
            1,
            2,
            3,
            8,
            15,
            48,
            105,
            384,
            945,
            3840,
            10395,
            46080,
            135135,
            645120,
            2027025,
        };

        assert(x <= 33 && "[Error]: double factorial overflows");
        if (x < Table.size())
            return Table[x];

        const uint64_t start = (x % 2 == 0) ? 14 : 15;
        uint64_t result = Table[start];
        for (uint64_t i = start + 2; i <= x; i += 2)
            result *= i;
        return result;
    }
}
