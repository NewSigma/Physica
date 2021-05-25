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
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"

using namespace Physica::Core;

void testLnGamma() {
    constexpr static int count = 2;
    constexpr static double value[count]{13.7, 0.3};
    constexpr static float floatResult[count]{21.77465, 1.095798};
    constexpr static double doubleResult[count]{21.77464517303463, 1.09579799481807552};
    
    for (int i = 0; i < count; ++i) {
        Scalar<Float, false> s(value[i]);
        auto temp = lnGamma(s);
        if ((fabs(float(temp) - floatResult[i]) / floatResult[i] >= 1E-6F))
            exit(EXIT_FAILURE);
    }

    for (int i = 0; i < count; ++i) {
        Scalar<Double, false> s(value[i]);
        auto temp = lnGamma(s);
        if ((fabs(double(temp) - doubleResult[i]) / doubleResult[i] >= 1E-14))
            exit(EXIT_FAILURE);
    }
}

void testGammaPQ() {
    using T = Scalar<Double, false>;
    constexpr static int count = 2;
    constexpr static double a[count]{13.7, 0.3};
    constexpr static double x[count]{2, 8};
    constexpr static double result[count]{5.309424005280372E-8, 0.99997576072630326};

    for (int i = 0; i < count; ++i) {
        auto temp = gammaP(T(a[i]), T(x[i]));
        if ((fabs(double(temp) - result[i]) >= result[i] * 1E-14))
            exit(EXIT_FAILURE);
        temp = gammaQ(T(a[i]), T(x[i]));
        if ((fabs(1 - double(temp) - result[i]) >= result[i] * 1E-9))
            exit(EXIT_FAILURE);
    }
}

void testBesselJ() {
    using T = Scalar<Double, false>;
    constexpr static int count = 2;
    constexpr static int n[count]{2, 5};
    constexpr static double x[count]{3, 3};
    constexpr static double result[count]{0.48609126058589107691, 0.043028434877047583925};

    for (int i = 0; i < count; ++i) {
        auto temp = besselJn(n[i], T(x[i]));
        if ((fabs(double(temp) - result[i]) >= result[i] * 1E-8))
            exit(EXIT_FAILURE);
    }
}

int main() {
    testLnGamma();
    testGammaPQ();
    testBesselJ();
    return 0;
}