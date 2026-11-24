/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Scalar/Complex.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void simd() {
        /* Div */ {
            auto x = SIMD<float64, 2>::random_uniform<Random<>>();
            auto y = SIMD<cfloat64, 2>::random_uniform<Random<>>();
            auto z = (x / y) * y;
            for (int i = 0; i < x.size(); ++i)
                expect<RandomSource>(scalarNear(x[i], z[i].real(), 2UL));
        }
        /* Tanh */ {
            auto x = SIMD<cfloat64, 2>::random_uniform<Random<>>();
            auto y = tanh(x);
            for (int i = 0; i < x.size(); ++i)
                expect<RandomSource>(tanh(x[i]) == y[i]);
        }
    }
}

int main() {
    static_assert(std::formattable<Complex<float64>, char>);
    simd();
    expect(softplus(cfloat64(-1000, 0)).isZero()); // Test underflow
    return 0;
}
