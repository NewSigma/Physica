/*
 * Copyright 2021-2024 WeiBo He.
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

#include <cassert>
#include <climits>
#include <cstdint>
#include <limits>
#include <algorithm>

namespace Physica::Core {
    /**
     * Calculate the number of arrangement $A_m^n$.
     * Using the definition, may be possible to optimize.
     */
    template<class IntType>
    IntType arrangement(IntType m, IntType n) {
        static_assert(std::numeric_limits<IntType>::is_integer, "[Error]: Invalid IntType");
        assert(m >= n && "[Error]: Invalid param");
        const IntType critical = m - n;
        IntType temp(m);
        IntType result = 1;
        while(temp > critical) {
            result *= temp;
            --temp;
        }
        return result;
    }
    /**
     * Calculate the number of combination $C_m^n$.
     * Using the definition, may be possible to optimize.
     */
    template<class IntType>
    IntType combination(IntType m, IntType n) {
        static_assert(std::numeric_limits<IntType>::is_integer, "[Error]: Invalid IntType");
        constexpr static int MaxNumberM = 16;
        constexpr static int K = (MaxNumberM - 4) / 2 + 1;
        constexpr static int TableLength = (MaxNumberM % 2 == 0) ? K * K : (K + 1) * K;
        constexpr static int64_t Table[TableLength]{
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
            120, 560, 1820, 4368, 8008, 11440, 12870
        };

        assert(m >= n && "[Error]: Invalid param");
        n = std::min(n, m - n);
        if (n == 0)
            return 1;
        if (n == 1)
            return m;
        const IntType k = (m - 4) / 2 + 1;
        const IntType rowShift = (m % 2 == 0) ? k * (k - 1) : k * k;
        return Table[rowShift + n - 2];
    }
}