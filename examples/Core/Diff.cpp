/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Scalar/Diff.h"

using namespace Physica;

static float64 forward() noexcept {
    using dfloat = Diff<float64, DiffMode::Forward>;
    dfloat x(0.5, 1); // x = 0.5, dx/dx = 1
    dfloat y = sin(x) / x; // Compute value as well as grad
    return y.grad();
}

static float64 reverse() noexcept {
    using dfloat = Diff<float64, DiffMode::Reverse>;
    dfloat x(0.5); // x = 0.5, initial grad is 0
    {
        auto y = sin(x) / x; // Compute value, sequence compute graph
        static_assert(std::same_as<decltype(y), CoDiff<dfloat>>); // return type contains a dfloat and a compute graph
        y.reverse(1); // We will propagate 1 reversely
    } // Start reverse propagate and destroy y
    return x.grad(); // x received grad
}

int main() {
    std::cout << std::format("{} {}\n", forward(), reverse());
    return 0;
}
