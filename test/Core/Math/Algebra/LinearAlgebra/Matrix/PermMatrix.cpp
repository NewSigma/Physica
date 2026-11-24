/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    void assign() {
        PermMatrix<T> perm({2, 1, 3, 0});
        perm = IdentityMatrix<T>(perm.getOrder());
        for (size_t i = 0; i < perm.getOrder(); ++i)
            expect(perm.getIndices()[i] == i);
    }

    void mulAndInv() {
        const PermMatrix<T> perm({0, 3, 2, 1});
        const VectorND<T> v{0, 1, 2, 3};
        const VectorND<T> perm_v = perm * v;
        {
            const VectorND<T> answer{0, 3, 2, 1};
            expect(vectorNear(perm_v, answer, 1E-15));
        }
        const auto inv = perm.inv();
        const VectorND<T> v1 = inv * perm_v;
        expect(vectorNear(v, v1, 1E-15));
    }

    void fromLAPACK() {
        PermMatrix<T> answer({2, 1, 3, 0});
        auto result = PermMatrix<T>::fromLAPACK(answer.toLAPACK<int64_t>());
        for (size_t i = 0; i < answer.getRow(); ++i)
            expect(result.getIndices()[i] == answer.getIndices()[i]);
    }
}

int main() {
    assign();
    mulAndInv();
    fromLAPACK();
    return 0;
}
