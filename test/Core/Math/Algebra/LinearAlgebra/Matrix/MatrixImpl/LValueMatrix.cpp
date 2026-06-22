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
    void diag() {
        // Test that diag() works for rectangular matrix
        using T = float32;
        using Matrix4x3 = DenseMatrix<T, MatrixMajor::Col, 4, 3>;
        auto compact = Matrix4x3::random_uniform<Random<>>(4, 3);
        auto corner = compact.topLeftCorner(3);
        expect(compact.diag().getLength() == 3);
        for (int i = 0; i < 3; ++i)
            expect(corner.diag()[i] == compact.diag()[i]);
    }

    void minorDiag() {
        // Test that 0 sub-diagonal is diagonal
        using T = float32;
        using Matrix4D = DenseMatrix<T, MatrixMajor::Col, 4, 4>;
        auto compact = Matrix4D::random_uniform<Random<>>(4, 4);
        auto corner = compact.topLeftCorner(3);
        for (int i = 0; i < 3; ++i)
            expect(corner.diag(0)[i] == compact.diag()[i]);
    }
}

int main() {
    diag();
    minorDiag();
    return 0;
}
