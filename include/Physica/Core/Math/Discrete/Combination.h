/*
 * Copyright 2021-2026 Weibo He.
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

#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <numeric>
#include "Physica/Core/Math/Calculus/SpecialFunctions.h"

namespace Physica {
    /**
     * Arrangement $A_m^n$.
     */
    [[nodiscard]] constexpr uint64_t arrangement(uint64_t m, uint64_t n) noexcept {
        assert(m >= n && "[Error]: Invalid param");
        uint64_t result = 1;
        for (uint64_t i = 0; i < n; ++i) {
            const uint64_t factor = m - i;
            assert(result <= std::numeric_limits<uint64_t>::max() / factor && "[Error]: Arrangement overflow");
            result *= factor;
        }
        return result;
    }
    /**
     * Combination $C_m^n$.
     */
    [[nodiscard]] constexpr uint64_t binomial(uint64_t m, uint64_t n) noexcept {
        [[maybe_unused]] constexpr static uint64_t MaxBinomialM = 67;
        constexpr static int MaxNumberM = 16;
        constexpr static int K = (MaxNumberM - 4) / 2 + 1;
        constexpr static int TableLength = (MaxNumberM % 2 == 0) ? K * K : (K + 1) * K;
        constexpr static std::array<int64_t, TableLength> Table = {
            6, // From 4
            10,
            15, 20,
            21, 35,
            28, 56, 70,
            36, 84, 126,
            45, 120, 210, 252,
            55, 165, 330, 462,
            66, 220, 495, 792, 924,
            78, 286, 715, 1287, 1716,
            91, 364, 1001, 2002, 3003, 3432,
            105, 455, 1365, 3003, 5005, 6435,
            120, 560, 1820, 4368, 8008, 11440, 12870,
        };

        assert(n <= m && "[Error]: Invalid param");
        assert(m <= MaxBinomialM && "[Error]: Binomial overflow");
        n = std::min(n, m - n);
        if (n == 0)
            return 1;
        if (n == 1)
            return m;

        if (m <= MaxNumberM) {
            const uint64_t k = (m - 4) / 2 + 1;
            const uint64_t rowShift = (m % 2 == 0) ? k * (k - 1) : k * k;
            return Table[rowShift + n - 2];
        }

        uint64_t result = 1;
        for (uint64_t i = 1; i <= n; ++i) {
            uint64_t numerator = m - n + i;
            uint64_t denominator = i;

            uint64_t g = std::gcd(numerator, denominator);
            numerator /= g;
            denominator /= g;

            g = std::gcd(result, denominator);
            result /= g;
            denominator /= g;

            assert(denominator == 1 && "[Error]: Binomial overflow");
            assert(result <= std::numeric_limits<uint64_t>::max() / numerator && "[Error]: Binomial overflow");
            result *= numerator;
        }
        return result;
    }

    template<Scalar T>
    T lnBinomial(T m, T n) noexcept {
        return lnGamma(m + T(1)) - lnGamma(n + T(1)) - lnGamma(m - n + T(1));
    }
}
