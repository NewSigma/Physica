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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;

namespace {
    void minorDiag() {
        // Test that 0 sub-diagonal is diagonal
        using T = float32;
        using Matrix4D = DenseMatrix<T, MatrixMajor::Col, 4, 4>;
        auto compact = Matrix4D::random_uniform<Random<>>(4, 4);
        auto corner = compact.topLeftCorner(3);
        for (int i = 0; i < 3; ++i)
            expect(abs_elem(corner).diag(0).calc(i) == abs_elem(compact).diag().calc(i));
    }

    void argmin() {
        const auto x = MatrixND<float32>::random_uniform<Random<>>(4, 5);
        auto index = x.argmin();
        expect(x.min() == x[index[0], index[1]]);
    }

    void sumRowCol() {
        const auto x = MatrixND<float32>::random_uniform<Random<>>(4, 5);
        VectorND<float32> buffer = x.sum_cols();
        for (size_t i = 0; i < x.getRow(); ++i)
            expect(buffer[i] == x.row(i).sum());

        buffer = x.sum_rows();
        for (size_t i = 0; i < x.getCol(); ++i)
            expect(buffer[i] == x.col(i).sum());
    }
}

static_assert(std::same_as<decltype(MatrixND<float32>{}.transpose().transpose()), MatrixND<float32>&&>, "[Error]: transpose-transpose should yield it self");

int main() {
    argmin();
    minorDiag();
    sumRowCol();
    return 0;
}
