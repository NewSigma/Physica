/*
 * Copyright 2021-2025 Weibo He.
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
#include "Functions.h"
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
