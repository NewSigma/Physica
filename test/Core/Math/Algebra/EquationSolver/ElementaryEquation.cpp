/*
 * Copyright 2021 WeiBo He.
 *
 * This file is part of Physica.

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
#include "Physica/Core/Math/Algebra/EquationSolver/ElementaryEquation.h"

int main() {
    using namespace Physica::Core;
    using ScalarType = Scalar<Double, false>;
    auto func = [](const ScalarType& x) { return exp(x) - ScalarType(5); };
    ScalarType result = bisectionMethod(func, ScalarType::Zero(), ScalarType(3));
    const ScalarType answer = ScalarType(1.6094379124341003746);
    if (!scalarNear(result, answer, 1E-15))
        return 1;
    return 0;
}
