/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/SparseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using SparseVec = SparseVector<T>;

static_assert(SparseVec::isSparse());
static_assert(!VectorND<T>::isSparse());

int main() {
    SparseVec v(4);
    v[0] = T(3);
    v[2] = T(7);
    expect(v.getNumNonZero() == 2);
    expect(v.getLength() == 4);
    expect(v.calc(0) == T(3));
    expect(v.calc(2) == T(7));
    expect(v.calc(1) == T(0));
    expect(v.calc(3) == T(0));
    return 0;
}
