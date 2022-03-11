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
#include "Physica/Utils/TestHelper.h"
#include "Functions.h"
#include "Physica/Core/Math/Optimization/Adam.h"

using namespace Physica::Utils;

using T = Scalar<Double, false>;

int main() {
    {
        Adam<T, Vector<T>> adam(Array<T, 6>({0.001, 0.9, 0.999, 1 - 1E-8, 1E-8, 1E-8}));
        adam.compute(func1<T>, Vector<T>{-1, -2, -5}, 0);
        if (!scalarNear(func1<T>(adam.getParams()), T(0), 1E-8))
            return 1;
    }
    return 0;
}