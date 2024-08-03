/*
 * Copyright 2023 Weibo He.
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

template<class T>
T func(T x, T y) {
    return square(x - T(1.0)) + square(y - T(2.0));
}

int main() {
    using PlainScalar = Scalar<Double>;
    {
        using ScalarType = Differentiable<PlainScalar, DiffMode::Forward, 1>;
        const PlainScalar x = 3;
        const PlainScalar y = 4;
        const ScalarType result = func(ScalarType(x, 1), ScalarType(y, 1));
        const PlainScalar answer = (x + y - 3.0) * 2.0;
        if (!scalarNear(result.getGrad(), answer, 1E-15))
            return 1;
    }
    {
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse, 1>;
        ScalarType x(3);
        ScalarType y(4);
        ScalarType result = func(x, y);
        result.reverse();
        if (!scalarNear(x.getGrad(), (x.getValue() - 1.0) * 2.0, 1E-15))
            return 1;
        if (!scalarNear(y.getGrad(), (y.getValue() - 2.0) * 2.0, 1E-15))
            return 1;
    }
    {
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse, 2>;
        ScalarType x(3);
        ScalarType result = square(x);
        result.reverse();
        if (!scalarNear(PlainScalar(x.getGrad<1>()), PlainScalar(6), 1E-15))
            return 1;
        x.getGrad().reverse();
        if (!scalarNear(x.getGrad<2>(), PlainScalar(2), 1E-15))
            return 1;
    }
    {
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse, 2>;
        ScalarType x(3);
        ScalarType y(4);
        ScalarType result = func(x, y);
        result.reverse();
        if (!scalarNear(PlainScalar(x.getGrad<1>()), PlainScalar(4), 1E-15))
            return 1;
        if (!scalarNear(PlainScalar(y.getGrad<1>()), PlainScalar(4), 1E-15))
            return 1;
        x.getGrad().reverse();
        if (!scalarNear(PlainScalar(x.getGrad<2>()), PlainScalar(2), 1E-15))
            return 1;
        if (!scalarNear(PlainScalar(y.getGrad<2>()), PlainScalar(0), 1E-15))
            return 1;
    }
    return 0;
}
