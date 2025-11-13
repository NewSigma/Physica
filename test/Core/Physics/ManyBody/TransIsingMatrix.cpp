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
#include "Physica/Core/Physics/ManyBody/Hamilton/TransIsingMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h"

using namespace Physica;
using T = float32;

namespace {
    template<BoundaryCond BC>
    void test() {
        TransIsingMatrix<T, 1, 3, BC> hamilton(1, 2, SquareLattice<1, BC>{{3}, 1});

        const auto size = hamilton.getRow();
        DenseMatrix<T> mat(size, size);
        for (size_t i = 0; i < size; ++i)
            mat.col(i) = hamilton * UnitVector<T>(i, size);

        if (!matrixNear(hamilton, mat, 1E-7))
            exit(EXIT_FAILURE);
    }
}

int main() {
    using enum BoundaryCond;
    test<OBC>();
    test<PBC>();
    return 0;
}
