/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using RandomSource = Random<>;

int main() {
    VectorND<float64> a(8, 1.0); // A vector of length 8, with its elements being 1.0
    VectorND<float64> b{1, 2, 3, 4, 5, 6, 7, 8}; // A vector of length 8, with its elements being 1..8
    auto c_ = VectorND<float32>::random_uniform<RandomSource>(8); // A random vector of length 8, with elements uniformly distributed on the interval (0, 1)
    VectorND<float64> c = c_; // Convert float32 vector to float64 vector
    VectorND<float64> d(8); // A vector of length 8, with its elements being undefined

    a[0] = 1; // Take element and assign value
    assert(a.getLength() == 8); // Length of vector is 8

    d = hadamard(sin(a) + b, square(c)); // Feel free to assemble your math operations
    std::cout << d << '\n'; // Output result
    return 0;
}
