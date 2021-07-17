/*
 * Copyright 2021 WeiBo He.
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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"

using namespace Physica::Core;

int main() {
    using T = Scalar<Double, false>;
    const Vector<T> x{2, 3, 4, 5};

    const size_t rank = x.getLength();
    Vector<T> v(rank);
    const T tau = x.houseHolder(v);
    const T beta = x[0].isNegative() ? v[0] : -v[0];
    v[0] = 1;

    Vector<T> result = x - tau * (v * x) * v;
    if (abs(T(result[0] - beta) / beta) > T(1E-15))
        return 1;
    for (size_t i = 1; i < rank; ++i)
        if (abs(result[i]) > T(1E-15))
            return 1;
    return 0;
}
