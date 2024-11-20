/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/MultiPrecision/Real.h"
#include "Physica/Core/MultiPrecision/Diff.h"

using namespace Physica::Core;

template<Scalar T>
T func(T x, T y) {
    return square(x - T(1.0)) + square(y - T(2.0));
}

using T = float64;

void testFunc() {
    bool good = true;
    {
        using ScalarType = Diff<T, DiffMode::Forward, 1>;
        const T x = 3;
        const T y = 4;
        const ScalarType result = func(ScalarType(x, 1), ScalarType(y, 1));
        const T answer = (x + y - 3.0) * 2.0;
        good &= scalarNear(result.getGrad(), answer, 1E-15);
    }
    {
        using ScalarType = Diff<T, DiffMode::Forward, 2>;
        ScalarType x{3, 1};
        ScalarType y = square(x);
        good &= scalarNear(y.template getGrad<2>(), float64(2), 1E-15);
    }
    {
        using ScalarType = Diff<T, DiffMode::Reverse, 1>;
        ScalarType x(3);
        ScalarType y(4);
        ScalarType result = func(x, y);
        result.reverse();
        good &= scalarNear(x.getGrad(), (x.getValue() - 1.0) * 2.0, 1E-15);
        good &= scalarNear(y.getGrad(), (y.getValue() - 2.0) * 2.0, 1E-15);
    }
    {
        using ScalarType = Diff<T, DiffMode::Reverse, 2>;
        ScalarType x(3);
        ScalarType result = square(x);
        result.reverse();
        good &= scalarNear(T(x.getGrad<1>()), T(6), 1E-15);
        x.getGrad().reverse();
        good &= scalarNear(x.getGrad<2>(), T(2), 1E-15);
    }
    {
        using ScalarType = Diff<T, DiffMode::Reverse, 2>;
        ScalarType x(3);
        ScalarType y(4);
        ScalarType result = func(x, y);
        result.reverse();
        good &= scalarNear(T(x.getGrad<1>()), T(4), 1E-15);
        good &= scalarNear(T(y.getGrad<1>()), T(4), 1E-15);
        x.getGrad().reverse();
        good &= scalarNear(T(x.getGrad<2>()), T(2), 1E-15);
        good &= scalarNear(T(y.getGrad<2>()), T(0), 1E-15);
    }
    if (!good)
        exit(EXIT_FAILURE);
}

void testMath() {
    using dfloat = Diff<T, DiffMode::Forward, 2>;
    bool good = true;
    dfloat x(3, 1);
    auto y = reciprocal(x);
    good &= scalarNear(y.getGrad().getValue(), -square(reciprocal(x.getValue())), 1E-15);
    good &= scalarNear(y.getGrad<2>(), pow(reciprocal(x.getValue()), T(3)) * T(2), 1E-15);

    y = sqrt(x);
    good &= scalarNear(y.getGrad().getValue(), reciprocal(T(2) * sqrt(x.getValue())), 1E-15);
    good &= scalarNear(y.getGrad<2>(), -reciprocal(T(4) * x.getValue() * sqrt(x.getValue())), 1E-15);
    if (!good)
        exit(EXIT_FAILURE);
}

void testSIMD() {
    T value[4]{1.5, -1.5, 0, 2};
    T grad1[4]{1, 1, 1, 0};
    T grad2[4]{0, 0, 0, 0};

    using dfloat = Diff<T, DiffMode::Forward, 2>;
    SIMD<dfloat, 4> packet{};
    packet.load({value, {grad1, grad2}});

    bool good = true;
    auto result = abs(packet);
    for (int i = 0; i < 4; ++i)
        good &= scalarNear(abs(dfloat(value[i], grad1[i])), result[i], 1E-15);

    result = square(packet);
    for (int i = 0; i < 4; ++i)
        good &= scalarNear(square(dfloat(value[i], grad1[i])), result[i], 1E-15);

    result = reciprocal(packet);
    for (int i = 0; i < 4; ++i) {
        if (value[i].isZero())
            continue;
        good &= scalarNear(reciprocal(dfloat(value[i], grad1[i])), result[i], 1E-15);
    }

    result = exp(packet);
    for (int i = 0; i < 4; ++i)
        good &= scalarNear(exp(dfloat(value[i], grad1[i])), result[i], 1E-15);

    if (!good)
        exit(EXIT_FAILURE);
}

int main() {
    testFunc();
    testMath();
    testSIMD();
    return 0;
}
