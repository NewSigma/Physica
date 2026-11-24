/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Calculus/Chebyshev.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Test.h"

using namespace Physica;
using T = float64;

namespace {
    void testFit(std::invocable<T> auto fn) {
        const int sample = 30;
        const T from(-1);
        const T to(1);

        VectorND<T> coeff(30);
        chebyshev_fit(from, to, coeff, fn);
        const T delta = (to - from) / T(sample);
        T x = from;
        for (int i = 0; i < sample; ++i) {
            x += delta;
            T answer = fn(x);
            T result = chebyshev_calc(from, to, coeff, x);
            expect(scalarNear(result, answer, 1E-14));
        }
    }

    void testFitEven(std::invocable<T> auto fn) {
        const int sample = 30;
        const T from(-1);
        const T to(1);

        VectorND<T> coeff(30);
        chebyshev_fit_even(from, to, coeff, fn);
        const T delta = (to - from) / T(sample);
        T x = from;
        for (int i = 0; i < sample; ++i) {
            x += delta;
            T answer = fn(x);
            T result = chebyshev_calc_even(from, to, coeff, x);
            expect(fabs(result.toMachine() - answer.toMachine()) < 1E-14);
        }
    }

    void testFitOdd(std::invocable<T> auto fn) {
        const int sample = 30;
        const T from(-1);
        const T to(1);

        VectorND<T> coeff(30);
        chebyshev_fit_odd(from, to, coeff, fn);
        const T delta = (to - from) / T(sample);
        T x = from;
        for (int i = 0; i < sample; ++i) {
            x += delta;
            T answer = fn(x);
            T result = chebyshev_calc_odd(from, to, coeff, x);
            expect(scalarNear(result, answer, 1E-14));
        }
    }
}

int main() {
    testFit([](const T& x) { return exp(x); });
    testFitEven([](const T& x) { return x * x; });
    testFitOdd([](const T& x) { return -x; });
    return 0;
}
