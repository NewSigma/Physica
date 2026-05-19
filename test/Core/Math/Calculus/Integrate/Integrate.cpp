/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using T = float64;

int main() {
    {
        auto simpleFn = [](T x) static noexcept -> T {
            return T(M_PI_2) * x * sin(T(M_PI) * x);
        };

        IntegrateRange<T, 1> range({-1}, {1});
        Integrate<IntegrateMethod::Rectangular, T, 1> rec(range, 0.01);
        expect(scalarNear(T(1), rec.solve(simpleFn), 1E-4));

        Integrate<IntegrateMethod::Ladder, T, 1> ladder(range, 0.01);
        expect(scalarNear(T(1), ladder.solve(simpleFn), 1E-4));

        Integrate<IntegrateMethod::Simpson, T, 1> simpson(range, 0.01);
        expect(scalarNear(T(1), simpson.solve(simpleFn), 1E-8));
    }
    {
        IntegrateRange<T, 1> range({1E-10}, {1});
        Integrate<IntegrateMethod::Tanh_Sinh, T, 1> tanh_sinh(range, 0.001, 3500);
        expect(scalarNear(T(23.025850929940456840), tanh_sinh.solve([](T x) -> T { return reciprocal(x); }), 1E-7));
    }
    {
        IntegrateRange<T, 2> range({0, 0}, {1, 1});
        Integrate<IntegrateMethod::MonteCarlo, T, 2> mc(range, 1000);

        const T answer = 1.317363136305819;
        using RandomSource = Random<MT19937, std::mt19937::default_seed>;
        T result, deviation;
        result = mc.solve_e<RandomSource>(6, [](Vector2D<T> x) -> T {
            return reciprocal(sqrt(square(x[0]) + sin(x[1])));
        }, deviation);
        expect(answer.isFinite() && deviation.isFinite() && (abs(answer - result) < deviation));
    }
    return 0;
}
