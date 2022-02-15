/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Utils/TestHelper.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;

int main() {
    const Vector<ScalarType> x{0, 1, 2};
    const Vector<ScalarType> y{5, -3, 2};
    for (size_t i = 0; i < x.getLength(); ++i)
        if (!scalarNear(lagrange(x, y, i), y[i], 1E-16))
            return 1;
    return 0;
}