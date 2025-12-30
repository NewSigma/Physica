/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Physics/ManyBody/Hamilton/AndersonMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"
#include "Test.h"

using namespace Physica;
using T = float32;

namespace {
    void test() {
        AndersonMatrix<T, 4, false> hamilton(1, 0, {0.2, 0.3, 0.4}, {0.2, 0.3, 0.4}, {3, 2});
        expect(hamilton.isSymm());

        const auto size = hamilton.getRow();
        DenseMatrix<T> mat(size, size);
        for (size_t i = 0; i < size; ++i)
            mat.col(i) = hamilton * UnitVector<T>(i, size);
        expect(matrixNear(hamilton, mat, 1E-7));
    }
}

int main() {
    test();
    return 0;
}
