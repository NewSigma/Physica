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
#include "Physica/Core/Math/Calculus/Chebyshev.h"
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

using namespace Physica::Core;
using T = Scalar<Double, false>;

template<class Function>
void testFit(Function func) {
    const int sample = 30;
    const T from(-1);
    const T to(1);

    Vector<T> coeff(30);
    chebyshev_fit(from, to, coeff, func);
    const T delta = (to - from) / T(sample);
    T x = from;
    for (int i = 0; i < sample; ++i) {
        x += delta;
        T answer = func(x);
        T result = chebyshev_calc(from, to, coeff, x);
        if (fabs(result.getTrivial() - answer.getTrivial()) > 1E-14) {
            std::cout << "testFit failed: " << result << ' ' << answer << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

template<class Function>
void testFitEven(Function func) {
    const int sample = 30;
    const T from(-1);
    const T to(1);

    Vector<T> coeff(30);
    chebyshev_fit_even(from, to, coeff, func);
    const T delta = (to - from) / T(sample);
    T x = from;
    for (int i = 0; i < sample; ++i) {
        x += delta;
        T answer = func(x);
        T result = chebyshev_calc_even(from, to, coeff, x);
        if (fabs(result.getTrivial() - answer.getTrivial()) > 1E-14) {
            std::cout << "testFitEven failed: " << x << result << ' ' << answer << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

template<class Function>
void testFitOdd(Function func) {
    const int sample = 30;
    const T from(-1);
    const T to(1);

    Vector<T> coeff(30);
    chebyshev_fit_odd(from, to, coeff, func);
    const T delta = (to - from) / T(sample);
    T x = from;
    for (int i = 0; i < sample; ++i) {
        x += delta;
        T answer = func(x);
        T result = chebyshev_calc_odd(from, to, coeff, x);
        if (fabs(result.getTrivial() - answer.getTrivial()) > 1E-14) {
            std::cout << "testFitOdd failed: " << result << ' ' << answer << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

int main() {
    testFit([](const T& x) { return exp(x); });
    testFitEven([](const T& x) { return x * x; });
    testFitOdd([](const T& x) { return -x; });
    return 0;
}
