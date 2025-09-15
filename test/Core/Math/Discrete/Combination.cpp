/*
 * Copyright 2024-2025 Weibo He.
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

using namespace Physica;

namespace {
    int64_t combination1(int m, int n) {
        if (n == 0 || m == n)
            return 1;
        return combination1(m - 1, n) + combination1(m - 1, n - 1);
    }
}

int main() {
    for (int m = 0; m <= 16; ++m)
        for (int n = 0; n <= m; ++n)
            if (combination1(m, n) != combination(m, n))
                return 1;
    return 0;
}
