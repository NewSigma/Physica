/*
 * Copyright 2021-2022 Weibo He.
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
#include "Functions.h"
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Math/Optimization/ConjugateGradient.h"
#include "Physica/Core/Math/Calculus/Differential.h"

using namespace Physica;

using ScalarType = float64;
using VectorType = Vector3D<ScalarType>;

template<class Functor>
VectorType grad(Functor func, const VectorType& at, ScalarType diffStep) {
    assert(diffStep.isPositive());
    const size_t length = at.getLength();
    VectorType result(length);
    VectorType copy = at;
    for (size_t i = 0; i < length; ++i) {
        const ScalarType buffer = copy[i];
        result[i] = Differential<ScalarType>::doublePoint([&](ScalarType alpha) {
            copy[i] = alpha;
            return func(copy);
        }, buffer, diffStep);
        copy[i] = buffer;
    }
    return result;
}

int main() {
    {
        auto func = func1<ScalarType>;
        ConjugateGradient<ScalarType, 3> cg(1);
        const ScalarType result = cg.solve(10, {-1, -2, -5}, func, [=](VectorType x) { return grad(func, x, 1E-5); });
        if (!scalarNear(result, ScalarType(0), 1E-14))
            return 1;
    }
    {
        auto func = func2<ScalarType>;
        ConjugateGradient<ScalarType, 3> cg(1);
        const ScalarType result = cg.solve(1E-15, {1, 3, 2}, func, [=](VectorType x) { return grad(func, x, 1E-6); });
        if (!scalarNear(result, ScalarType(2.25), 1E-14))
            return 1;
    }
    {
        auto func = rosenbrock<ScalarType>;
        ConjugateGradient<ScalarType, 3> cg(1);
        const ScalarType result = cg.solve(1E-13, {-1.2, 1}, func, [=](VectorType x) { return grad(func, x, 3E-6); });
        if (!scalarNear(result, ScalarType(0), 1E-18))
            return 1;
    }
    return 0;
}
