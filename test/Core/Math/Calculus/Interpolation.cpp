/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

namespace {
    void testLaglange() {
        const VectorType x{0, 1, 2};
        const VectorType y{5, -3, 2};
        for (size_t i = 0; i < x.getLength(); ++i)
            expect(scalarNear(lagrange(x, y, i), y[i], 1E-16));
    }

    void testFFT1D() {
        const auto data = VectorType::template random_normal<RandomSource>(20);
        const auto result = interpolate_fft(data, 100);

        const size_t delta = result.getLength() / data.getLength();
        for (size_t i = 0; i < result.getLength(); i += delta)
            expect(scalarNear(data[i / delta], result[i], 1E-11));

        for (size_t i = 0; i < data.getLength(); ++i) {
            ScalarType result = interpolate_fft(data, i, data.getLength());
            expect(scalarNear(result, data[i], 1E-11));
        }
        /* Test periodicity */ {
            const ScalarType value = interpolate_fft(data, 0.2, 1);
            const ScalarType value1 = interpolate_fft(data, 1.2, 1);
            expect(scalarNear(value, value1, 1E-12));
        }
    }

    void testFFT3D() {
        using TensorType = DenseTensor<ScalarType, 3>;
        const auto data = TensorType::random_uniform<RandomSource>({5, 5, 5});
        {
            const auto result = interpolate_fft(data, {10, 10, 10});
            for (size_t i = 0; i < data.getSize(); ++i) {
                auto indices = data.toIndexND(i);
                auto indices1 = indices;
                for (auto& id : indices1)
                    id *= 2;
                expect(scalarNear(data.asArray()[i], result[indices1], 1E-8));
            }
        }
        /* Test periodicity */ {
            const ScalarType value = interpolate_fft(data, {0.2, 0.1, 0.3}, {1, 1, 1});
            const ScalarType value1 = interpolate_fft(data, {1.2, 0.1, 0.3}, {1, 1, 1});
            expect(scalarNear(value, value1, 1E-12));
        }
    }
}

int main() {
    testLaglange();
    testFFT1D();
    testFFT3D();
    return 0;
}
