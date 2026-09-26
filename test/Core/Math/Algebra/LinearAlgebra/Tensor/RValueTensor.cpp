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
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    void fiber() {
        const auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        const auto y = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        static_assert(!(x + y).isLValueTensor());

        auto fiber = (x + y).fiber(1, var(), 2);
        static_assert(!fiber.isLValueVector());
        expect(VectorND<T>(fiber) == fiber);
    }

    void slice() {
        const auto x = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        const auto y = Tensor3D<T>::random_uniform<Random<>>({4, 4, 4});
        auto slice = (x + y).slice(1, var(), var());
        static_assert(!slice.isLValueMatrix());
        expect(MatrixND<T>(slice) == slice);
    }
}

int main() {
    fiber();
    slice();
    return 0;
}
