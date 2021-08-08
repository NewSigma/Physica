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
#include "Functions.h"
#include "Physica/Core/Math/Optimization/ConjugateGradient.h"

using namespace Physica::Core::Math;

using T = Scalar<Double, false>;

int main() {
    {
        ConjugateGradient cg(func1<T>, Vector<T>{-1, -2, -5}, T(1E-14), T(0.00001));
        if (fabs(cg.compute().getTrivial()) > 1E-14)
            return 1;
    }
    {
        ConjugateGradient cg(func2<T>, Vector<T>{1, 3, 2}, T(1E-14), T(0.0001));
        if (fabs(cg.compute().getTrivial() - 2.25) / 2.25 > 1E-14)
            return 1;
    }
    return 0;
}
