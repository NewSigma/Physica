/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Test.h"

using namespace Physica;

int main() {
    auto x = DenseTensor<float64, 3>::random_uniform<Random<>>({4, 4, 4});
    auto fiber = x.fiber(1, {1, Dynamic, 2});
    for (int i = 0; i < fiber.getLength(); ++i)
        expect(x[1, i, 2] == fiber[i]);

    auto slice = x.slice(1, 2, {1, Dynamic, Dynamic});
    for (int r = 0; r < fiber.getLength(); ++r)
        for (int c = 0; c < fiber.getLength(); ++c)
        expect(x[1, r, c] == slice[r, c]);
    return 0;
}
