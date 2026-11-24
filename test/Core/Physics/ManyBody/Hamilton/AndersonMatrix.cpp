/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
