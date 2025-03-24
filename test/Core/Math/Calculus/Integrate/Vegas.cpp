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
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Physica/Core/Math/Calculus/Integrate/Vegas.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using RandomSource = Random<MT19937, 12345>;
using T = float64;
/**
 * Reference:
 * [1] J. Comput. Phys. 439, 110386 (2021); https://doi.org/10.1016/j.jcp.2021.110386
 */
int main() {
    const Vector4D<T> r1{1.0 / 3, 0.5, 0.5, 0.5};
    const Vector4D<T> r2{2.0 / 3, 0.5, 0.5, 0.5};
    /* Eq.26 of [1] */ {
        auto func = [&](const VectorND<T>& x) {
            return exp(T(-100) * (x - r1).squaredNorm()) + exp(T(-100) * (x - r2).squaredNorm());
        };
        auto vegas = Vegas<T, false>({0, 0, 0, 0}, {1, 1, 1, 1}, 50, 10000);
        vegas.integral<decltype(func), RandomSource>(func);

        const T temp = erf(T(5));
        const T answer = T(std::numbers::pi * std::numbers::pi / 10000) * square(temp) * temp * (erf(T(10.0 / 3)) + erf(T(20.0 / 3)));
        const T result = vegas.calcMean();
        if (abs(answer - result) > T(2) * vegas.calcDevia())
            return 1;
    }
    /* Eq.28 of [1] */ {
        constexpr double R = 2.0 / 30;
        auto func = [&](const VectorND<T>& x) {
            const bool flag1 = (x - r1).squaredNorm() < T(R * R);
            const bool flag2 = (x - r2).squaredNorm() < T(R * R);
            return T((flag1 || flag2) ? 1 : 0);
        };
        auto vegas = Vegas<T, false>({0, 0, 0, 0}, {1, 1, 1, 1}, 50, 100000, 1000, 0.2);
        vegas.integral<decltype(func), RandomSource>(func);

        const T answer = square(T(std::numbers::pi * R * R));
        const T result = vegas.calcMean();
        if (abs(answer - result) > T(2) * vegas.calcDevia())
            return 1;
    }
    {
        auto func = [&](const VectorND<T>& x) {
            return T(-100) * (x - r1).squaredNorm();
        };
        auto vegas = Vegas<T, true>({0, 0, 0, 0}, {1, 1, 1, 1}, 50, 10000);
        vegas.integral<decltype(func), RandomSource>(func);

        const T temp = erf(T(5));
        const T answer = T(std::numbers::pi * std::numbers::pi / 20000) * square(temp) * temp * (erf(T(10.0 / 3)) + erf(T(20.0 / 3)));
        const T result = exp(vegas.calcLnMean());
        if (abs(answer - result) > T(2) * exp(vegas.calcLnDevia()))
            return 1;
    }
    return 0;
}
