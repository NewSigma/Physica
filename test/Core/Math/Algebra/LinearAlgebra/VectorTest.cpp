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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;

int main() {
    using T = Scalar<Float, false>;
    Vector<T> v{1, 2, 3, 4};
    auto col = v.copyToColMatrix();
    auto row = v.moveToRowMatrix();

    Vector<Scalar<Double, false>> v1{3.845971,0.000000,0.000000};
    Vector<Scalar<Double, false>> v2{-0.007733,3.835502,0.000000};
    Vector<Scalar<Double, false>> v3(v1.crossProduct(v2));
    if (fabs(v3.norm().getTrivial() / 2 - 7.375614) / 7.375614 > 10E-8)
        return 1;
    return 0;
}
