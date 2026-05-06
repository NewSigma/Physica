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

namespace {
    void fiber(const DenseTensor<float64, 3>& x) {
        auto fiber = x.fiber(1, {1, Dynamic, 2});
        for (int i = 0; i < fiber.getLength(); ++i)
            expect(x[1, i, 2] == fiber[i]);

        VectorND<float64> v = fiber;
        expect(v == fiber);
    }

    void slice(const DenseTensor<float64, 3>& x) {
        auto slice = x.slice(1, 2, {1, Dynamic, Dynamic});
        for (int r = 0; r < x.dim(1); ++r)
            for (int c = 0; c < x.dim(2); ++c)
            expect(x[1, r, c] == slice[r, c]);

        MatrixND<float64> m = slice;
        expect(m == slice);
    }
}

int main() {
    auto x = DenseTensor<float64, 3>::random_uniform<Random<>>({4, 4, 4});
    fiber(x);
    slice(x);
    return 0;
}
