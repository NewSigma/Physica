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
#include <iostream>
#include "Physica/Core/Math/Transform/DFT.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;

ScalarType func(ScalarType x) {
    if (x < ScalarType::One())
        return ScalarType::One();
    else
        return -ScalarType::One();
}

int main() {
    const size_t N = 50;
    const Vector<ScalarType> data = Vector<ScalarType>::linspace(ScalarType::Zero(), ScalarType::Two(), N);
    DFT<ScalarType> dft(data, ScalarType(2.0 / N));
    dft.transform();
    dft.invTransform();
    Vector<ScalarType> data1 = toRealVector(dft.getData());
    std::cout << (data1) << std::endl;
    return 0;
}
