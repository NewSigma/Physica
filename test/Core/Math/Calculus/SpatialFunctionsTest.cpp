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
#include "Physica/Core/Math/Calculus/SpatialFunctions.h"

using namespace Physica::Core;

int main() {
    constexpr static int count = 2;
    constexpr static double value[count]{13.7, 0.3};
    constexpr static float floatResult[count]{21.77465, 1.095798};
    constexpr static double doubleResult[count]{21.77464517303463, 1.09579799481807552};
    
    for (int i = 0; i < count; ++i) {
        Scalar<Float, false> s(value[i]);
        auto temp = lnGamma(s);
        if ((fabs(float(double(temp)) - floatResult[i]) / floatResult[i] >= 1E-6F))
            return 1;
    }

    for (int i = 0; i < count; ++i) {
        Scalar<Double, false> s(value[i]);
        auto temp = lnGamma(s);
        if ((fabs(double(temp) - doubleResult[i]) / doubleResult[i] >= 1E-14))
            return 1;
    }
    return 0;
}