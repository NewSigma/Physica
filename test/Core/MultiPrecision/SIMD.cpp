/*
 * Copyright 2024 Weibo He.
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
#include <iostream>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/MultiPrecision/Complex.h"

using namespace Physica::Core;
using RandomType = Random<MT19937, 10000>;

int main() {
    const auto x = SIMD<cfloat32, 2>::random_uniform(RandomType::getInstance());
    const auto result = sqrt(x);
    for (int i = 0; i < 2; ++i) {
        if (!scalarNear(result[i], sqrt(x[i]), 1E-7))
            return 1;
    }
    return 0;
}
