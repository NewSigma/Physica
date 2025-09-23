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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using RandomSource = Random<MT19937>;
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
