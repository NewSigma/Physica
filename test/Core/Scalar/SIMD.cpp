/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Complex.h"

using namespace Physica;
using RandomSource = Random<MT19937, 10000>;

template<Scalar T, int Size>
void testMath(const SIMD<T, Size>& x) {
    constexpr double prec = T::Prec == Float32 ? 1E-6 : 1E-15;
    /* Divide */ {
        const auto y = SIMD<T, Size>::template random_uniform<RandomSource>();
        const auto result = x / y;
        for (int i = 0; i < Size; ++i) {
            if (!scalarNear(result[i], x[i] / y[i], prec))
                exit(EXIT_FAILURE);
        }
    }
    /* Exp */ {
        const auto result = exp(x);
        for (int i = 0; i < Size; ++i) {
            if (!scalarNear(result[i], exp(x[i]), prec))
                exit(EXIT_FAILURE);
        }
    }
    /* Ln */ {
        const auto result = ln(x);
        for (int i = 0; i < Size; ++i) {
            if (!scalarNear(result[i], ln(x[i]), prec))
                exit(EXIT_FAILURE);
        }
    }
    /* LnCosh */ {
        const auto result = lncosh(x);
        for (int i = 0; i < Size; ++i) {
            if (!scalarNear(result[i], lncosh(x[i]), prec))
                exit(EXIT_FAILURE);
        }
    }
}

template<Scalar T, int Size>
void test() {
    constexpr double prec = T::Prec == Float32 ? 1E-6 : 1E-15;
    auto x = SIMD<T, Size>(std::numeric_limits<T>::min());
    /* Divide */ {
        const auto result = SIMD<T, Size>(0) / x;
        for (int i = 0; i < Size; ++i) {
            if (!result[i].isZero())
                exit(EXIT_FAILURE);
        }
    }
    testMath<T, Size>(x);

    x = SIMD<T, Size>::template random_uniform<RandomSource>();
    testMath<T, Size>(x);
    /* Sqrt */ {
        const auto result = sqrt(x);
        for (int i = 0; i < Size; ++i) {
            if (!scalarNear(result[i], sqrt(x[i]), prec))
                exit(EXIT_FAILURE);
        }
    }
}

int main() {
    test<cfloat32, 2>();
    test<cfloat32, 4>();
    test<cfloat64, 1>();
    test<cfloat64, 2>();
    return 0;
}
