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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    void fiber() {
        auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        static_assert(x.fiber(1, 2, var()).isCompact());

        auto fiber = x.fiber(1, var(), 2);
        static_assert(!decltype(fiber)::isCompact());
        static_assert(decltype(fiber)::getStrideAtCompile() == Dynamic);
        expect(fiber.getStride() == x.getStride(1));
        for (int i = 0; i < fiber.getLength(); ++i)
            expect(x[1, i, 2] == fiber[i]);
    }

    void slice() {
        auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        auto slice = x.slice(1, var(), var());
        static_assert(slice.isCompact());
        static_assert(slice.getMajor() == MatrixMajor::Row);
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                expect(x[1, r, c] == slice[r, c]);

        MatrixND<T> m = slice;
        expect(m == slice);
    }
}

int main() {
    fiber();
    slice();
    return 0;
}
