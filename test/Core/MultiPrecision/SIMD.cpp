/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/MultiPrecision/Complex.h"

using namespace Physica::Core;
using RandomType = Random<MT19937, 10000>;

template<Scalar T, int Size>
void test() {
    constexpr double prec = T::Option == Float32 ? 1E-6 : 1E-15;
    const auto x = SIMD<T, Size>::template random_uniform<RandomType>();
    /* Divide */ {
        const auto y = SIMD<T, Size>::template random_uniform<RandomType>();
        const auto result = x / y;
        for (int i = 0; i < Size; ++i) {
            if (!scalarNear(result[i], x[i] / y[i], prec))
                exit(EXIT_FAILURE);
        }
    }
    /* Sqrt */ {
        const auto result = sqrt(x);
        for (int i = 0; i < Size; ++i) {
            if (!scalarNear(result[i], sqrt(x[i]), prec))
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

int main() {
    test<cfloat32, 2>();
    test<cfloat32, 4>();
    test<cfloat64, 1>();
    test<cfloat64, 2>();
    return 0;
}
