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
#include "Functions.h"
#include "Physica/Core/Math/Optimization/ConjugateGradient.h"

using namespace Physica::Core::Math;

using ScalarType = Scalar<Double, false>;

int main() {
    {
        ConjugateGradient cg(func1<ScalarType>, Vector<ScalarType>{-1, -2, -5}, ScalarType(10), ScalarType(1), ScalarType(1E-5));
        if (!scalarNear(cg.compute(), ScalarType::Zero(), 1E-15))
            return 1;
    }
    {
        ConjugateGradient cg(func2<ScalarType>, Vector<ScalarType>{1, 3, 2}, ScalarType(1E-15), ScalarType(1), ScalarType(1E-6));
        if (!scalarNear(cg.compute(), ScalarType(2.25), 1E-14))
            return 1;
    }
    {
        ConjugateGradient cg(rosenbrock<ScalarType>, Vector<ScalarType>{-1.2, 1}, ScalarType(1E-13), ScalarType(1), ScalarType(3E-6));
        if (!scalarNear(cg.compute(), ScalarType::Zero(), 1E-19))
            return 1;
    }
    return 0;
}
