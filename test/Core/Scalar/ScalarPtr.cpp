/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Scalar/Diff.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using dfloat = Diff<float64, DiffMode::Forward, 2>;

int main() {
    // Test ScalarPtr is a tuple
    auto [p1, p2, p3] = ScalarPtr<dfloat>{};
    static_assert(std::same_as<decltype(p1), T*>);
    static_assert(std::same_as<decltype(p3), T*>);
    return 0;
}
