/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Optimization/TestFunction.h"
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Math/Optimization/ConjugateGradient.h"
#include "Test.h"

using namespace Physica;

using T = float64;

namespace {
    template<Vector V>
    V grad(std::invocable<V> auto fn, const V& at, T diffStep) {
        assert(diffStep.isPositive());
        const size_t length = at.getLength();
        V result(length);
        V copy = at;
        for (size_t i = 0; i < length; ++i) {
            const T buffer = copy[i];
            result[i] = Differential<T>::doublePoint([&](T alpha) {
                copy[i] = alpha;
                return fn(copy);
            }, buffer, diffStep);
            copy[i] = buffer;
        }
        return result;
    }
}

int main() {
    {
        auto func = func1<T>;
        ConjugateGradient<T, 3> cg(1);
        const T result = cg.solve(10, {-1, -2, -5}, func, [=](Vector3D<T> x) { return grad(func, x, 1E-5); });
        expect(scalarNear(result, T(0), 1E-14));
    }
    {
        auto func = func2<T>;
        ConjugateGradient<T, 3> cg(1);
        const T result = cg.solve(1E-15, {1, 3, 2}, func, [=](Vector3D<T> x) { return grad(func, x, 1E-6); });
        expect(scalarNear(result, T(2.25), 1E-14));
    }
    {
        auto func = rosenbrock<T>;
        ConjugateGradient<T, 2> cg(1);
        const T result = cg.solve(1E-13, {-1.2, 1}, func, [=](Vector2D<T> x) { return grad(func, x, 3E-6); });
        expect(scalarNear(result, T(0), 1E-18));
    }
    return 0;
}
