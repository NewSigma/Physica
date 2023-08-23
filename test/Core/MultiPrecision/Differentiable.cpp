/*
 * Copyright 2023 WeiBo He.
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
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/MultiPrecision/Differentiable.h"

using namespace Physica::Core;
using ScalarType = Differentiable<Scalar<Double>>;

template<class T>
T func(T x, T y) {
    using PlainType = typename T::PlainScalar;
    return square(x - PlainType(1.0)) + square(y - PlainType(2.0));
}

int main() {
    const Scalar<Double> x = 3;
    const Scalar<Double> y = 4;
    const ScalarType result = func(ScalarType(x, 1), ScalarType(y, 1));
    const Scalar<Double> answer = (x + y - 3.0) * 2.0;
    if (!scalarNear(result.getTangent(), answer, 1E-15))
        return 1;
    return 0;
}
