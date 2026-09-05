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
    void forward() {
        using T = float64;
        using dfloat = Diff<T, DiffMode::Forward>;
        auto m = MatrixND<dfloat>::random_uniform<RandomSource>(4, 4);
        auto x = VectorND<dfloat>::random_uniform<RandomSource>(4);
        auto y = VectorND<dfloat>::junk(4);
        y = m * x;
        expect<RandomSource>(vectorNear(y.values(), VectorND<T>(m.values() * x.values()), 1E-13));
        expect(y.grads().isZero());
    }

    void reverse() {
        using ScalarType = Diff<float32, DiffMode::Reverse>;
        using MatrixType = DenseMatrix<float32, MatrixMajor::Col, 3, 3>;
        using DVector = Vector3D<ScalarType>;
        using DMatrix = DenseMatrix<ScalarType, MatrixMajor::Col, 3, 3>;
        DMatrix m{1, 2, 3, 4, 5, 6, 7, 8, 9};
        DVector x{1, 2, 3};
        {
            CoDiff<DVector> y = m * x;
            y.sum().reverse();
        }
        const Vector3D<float32> dx = x.grads();
        for (size_t i = 0; i < dx.getLength(); ++i)
            expect(scalarNear(dx[i], m.values().col(i).sum(), 1E-15));

        const MatrixType dm = m.grads();
        for (size_t r = 0; r < dm.getRow(); ++r)
            for (size_t c = 0; c < dm.getCol(); ++c)
                expect(scalarNear(dm[r, c], x.values().calc(c), 1E-15));
    }
}

int main() {
    forward();
    reverse();
    return 0;
}
