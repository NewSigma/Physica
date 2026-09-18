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
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using TensorType = Tensor3D<T>;
using RandomSource = Random<>;

namespace {
    void compact() {
        const TensorType x = TensorType::random_uniform<RandomSource>({4, 4, 4});
        expect(x.data() == x.asArray().data());
        expect(x.data() == x.data_handle());
        expect(x.data_ptr({1, 2, 3}) == x.data() + x.toIndex1D({1, 2, 3}));

        for (size_t i = 0; i < x.dim(0); ++i)
            for (size_t j = 0; j < x.dim(1); ++j)
                for (size_t k = 0; k < x.dim(2); ++k)
                    expect(x.data_ptr({i, j, k}) == &x[i, j, k]);
    }
}

static_assert(TensorType::isCompact(), "DenseTensor is a compact object");

int main() {
    compact();
    return 0;
}
