/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Math/Calculus/Interpolation.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;

void testLaglange() {
    const VectorType x{0, 1, 2};
    const VectorType y{5, -3, 2};
    for (size_t i = 0; i < x.getLength(); ++i)
        if (!scalarNear(lagrange(x, y, i), y[i], 1E-16))
            exit(EXIT_FAILURE);
}

void testFFT1D() {
    std::mt19937 gen{};
    const auto data = VectorType::random_normal(20, gen);
    const auto result = interpolate_fft(data, 100);

    const size_t delta = result.getLength() / data.getLength();
    for (size_t i = 0; i < result.getLength(); i += delta) {
        if (!scalarNear(data[i / delta], result[i], 1E-12))
            exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < data.getLength(); ++i) {
        ScalarType result = interpolate_fft(data, i, data.getLength());
        if (!scalarNear(result, data[i], 1E-13))
            exit(EXIT_FAILURE);
    }
}

void testFFT3D() {
    std::mt19937 gen{};
    const auto data = VectorType::random_normal(20, gen);
    const auto result = interpolate_fft(data, 100);
}

int main() {
    testLaglange();
    testFFT1D();
    testFFT3D();
}