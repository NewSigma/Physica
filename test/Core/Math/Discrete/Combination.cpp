/*
 * Copyright 2024-2026 Weibo He.
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
#include "Physica/Core/Math/Discrete/Combination.h"
#include "Test.h"

using namespace Physica;

static_assert(arrangement(0, 0) == 1);
static_assert(arrangement(5, 2) == 20);
static_assert(arrangement(20, 20) == 2432902008176640000ULL);

static_assert(binomial(5, 2) == 10);
static_assert(binomial(16, 8) == 12870);
static_assert(binomial(17, 8) == 24310);
static_assert(binomial(67, 33) == 14226520737620288370ULL);
static_assert(binomial(67, 66) == 67);

namespace {
    void combination() {
        auto recursive = [](this auto self, int m, int n) -> uint64_t {
            if (n == 0 || m == n)
                return 1;
            return self(m - 1, n) + self(m - 1, n - 1);
        };

        for (int m = 0; m <= 16; ++m)
            for (int n = 0; n <= m; ++n)
                expect(recursive(m, n) == binomial(m, n));
    }
}

int main() {
    combination();
    return 0;
}
