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
#include "TestHelper.h"
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"

int main() {
    using namespace Physica::Core;
    using ScalarType = Scalar<Double, false>;
    {
        auto func = [](const ScalarType& x) { return hermiteH(5, x); };
        auto result = Differential<ScalarType>::forward(func, ScalarType(3), ScalarType(1E-8));
        const ScalarType answer = ScalarType(8760);
        if (!scalarNear(result, answer, 1E-8))
            return 1;

        result = Differential<ScalarType>::backward(func, ScalarType(3), ScalarType(1E-8));
        if (!scalarNear(result, answer, 1E-7))
            return 1;

        result = Differential<ScalarType>::doublePoint(func, ScalarType(3), ScalarType(1E-5));
        std::cout << result << std::endl;
        if (!scalarNear(result, answer, 1E-10))
            return 1;
    }
    return 0;
}
