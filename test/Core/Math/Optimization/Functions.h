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
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

using namespace Physica::Core;

template<class T>
T func1(const Vector<T>& v) {
    const T& x = v[0];
    const T& y = v[1];
    const T& z = v[2];
    return x * x + y * y + z * z - x * y - x * z - y * z;
}

template<class T>
T func2(const Vector<T>& v) {
    const T& x = v[0];
    const T& y = v[1];
    const T& z = v[2];
    const T term1 = x + y;
    const T term2 = y + z;
    const T term3 = x + z;
    return (reciprocal(term1 * term1) + reciprocal(term2 * term2) + reciprocal(term3 * term3)) * T(x * y + x * z + y * z);
}
