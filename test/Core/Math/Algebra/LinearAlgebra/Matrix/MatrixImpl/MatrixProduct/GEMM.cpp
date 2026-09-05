/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Scalar/Diff.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void simple2D() {
        using M = DenseMatrix<float32, MatrixMajor::Col, 2, 2>;
        const M m1{{1, 2}, {2, 1}};
        const M m2{{3, 3}, {1, 5}};
        const M result = m1 * m2;
        const M answer{{9, 9}, {11, 7}};
        expect(matrixNear(result, answer, std::numeric_limits<float>::epsilon()));
    }

    void forward() {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Forward>;
        auto a = MatrixND<dfloat>::random_uniform<RandomSource>(4, 4);
        auto b = MatrixND<dfloat>::random_uniform<RandomSource>(4, 4);
        a.grads().random_uniform<RandomSource>();
        b.grads().random_uniform<RandomSource>();
        auto c = MatrixND<dfloat>::junk(4, 4);
        c = a * b;
        expect<RandomSource>(matrixNear(c.values(), MatrixND<T>(a.values() * b.values()), 1E-13));
    }
}

int main() {
    simple2D();
    forward();
    return 0;
}
