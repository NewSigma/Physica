/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

int main() {
    VectorND<float64> a(8, 1.0); // A vector of length 8, with its elements are 1.0
    VectorND<float64> b{1, 2, 3, 4, 5, 6, 7, 8}; // A vector of length 8, with its elements are 1..8
    auto c_ = VectorND<float32>::random_uniform<RandomSource>(8); // A random vector of length 8, with elements uniformly distributed on the interval (0, 1)
    VectorND<float64> c = c_; // Convert float32 vector to float64 vector
    VectorND<float64> d(8); // A vector of length 8, with its elements undefined

    a[0] = 1; // Take element and assign value
    assert(a.getLength() == 8); // Length of vector is 8

    d = hadamard(sin(a) + b, square(c)); // Feel free to assemble your math operations
    std::cout << d.format() << std::endl; // Output result
    return 0;
}
