/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/J1J2Matrix.h"
#include "Test.h"

using namespace Physica;
using T = float32;

namespace {
    template<BoundaryCond BC>
    void test() {
        J1J2Matrix<T, 2, 9, BC> hamilton(1, 0.3, SquareLattice<2, BC>{{3, 3}, 1});
        expect(hamilton.isSymm());

        const auto size = hamilton.getRow();
        DenseMatrix<T> mat(size, size);
        for (size_t i = 0; i < size; ++i)
            mat.col(i) = hamilton * UnitVector<T>(i, size);
        expect(matrixNear(hamilton, mat, 1E-7));
    }
}

int main() {
    using enum BoundaryCond;
    test<PBC>();
    return 0;
}
