/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using RandomSource = Random<>;
using DVector = device_obj<VectorND<float64>>; // Simply dress anything with device_obj<> to get its CUDA version

int main() {
    DVector a(8, 1.0);
    auto b_ = VectorND<float64>{1, 2, 3, 4, 5, 6, 7, 8}; // Any host vector
    DVector b = b_.toDeviceAsync(); // Async version, pass to GPU
    auto c = DVector::random_uniform<RandomSource>(8);
    DVector d(8);

    d = hadamard(sin(a) + b, square(c)); // For any combination, the operators are fused and only a single CUDA kernel is issued.
    std::cout << d.toHost().format() << '\n'; // Sync version, wait for the GPU to finish
    return 0;
}
