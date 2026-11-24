/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
