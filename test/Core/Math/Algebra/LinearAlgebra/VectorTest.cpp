/*
 * Copyright 2020-2021 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"

using namespace Physica::Core;

int main() {
    using T = Scalar<Float, false>;
    Vector<T> v1{3.845971,0.000000,0.000000};
    Vector<T> v2{-0.007733,3.835502,0.000000};
    Vector<T> v3(v1.crossProduct(v2));
    if (!scalarNear(v3.norm() / T::Two(), T(7.375614), 1E-7))
        return 1;
    return 0;
}
