/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;

ScalarType func(ScalarType x) {
    return ScalarType(M_PI_2) * x * sin(ScalarType(M_PI) * x);
}

int main() {
    {
        IntegrateRange<ScalarType, 1> range({-1}, {1});
        Integrate<Rectangular, ScalarType, 1> rec(range, 0.01);
        if (!scalarNear(ScalarType::One(), rec.solve(func), 1E-4))
            return 1;

        Integrate<Ladder, ScalarType, 1> ladder(range, 0.01);
        if (!scalarNear(ScalarType::One(), ladder.solve(func), 1E-4))
            return 1;

        Integrate<Simpson, ScalarType, 1> simpson(range, 0.01);
        if (!scalarNear(ScalarType::One(), simpson.solve(func), 1E-8))
            return 1;
    }
    {
        IntegrateRange<ScalarType, 1> range({1E-10}, {1});
        Integrate<Tanh_Sinh, ScalarType, 1> tanh_sinh(range, 0.001, 3500);
        if (!scalarNear(ScalarType(23.025850929940456840), tanh_sinh.solve([](ScalarType x) -> ScalarType { return reciprocal(x); }), 1E-7))
            return 1;
    }
    return 0;
}