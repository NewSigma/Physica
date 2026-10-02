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
