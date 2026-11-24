/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/TransIsingMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h"
#include "Test.h"

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
        expect(matrixNear(hamilton, mat, 1E-7));
    }
}

int main() {
    using enum BoundaryCond;
    test<OBC>();
    test<PBC>();
    return 0;
}
