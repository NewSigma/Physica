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
#include "TestHelper.h"
#include "Physica/Core/Physics/ElectronStructure/HP/GaussBase.h"

using namespace Physica::Core;
using namespace Physica::Core::Physics;

using ScalarType = Scalar<Double, false>;

ScalarType overlap_1s_1s(const ScalarType& alpha1,
                         const Vector<ScalarType, 3>& v1,
                         const ScalarType& alpha2,
                         const Vector<ScalarType, 3>& v2) {
    const ScalarType alpha_sum = alpha1 + alpha2;
    const ScalarType temp = ScalarType(M_PI) / alpha_sum;
    const ScalarType factor = temp * sqrt(temp);
    return factor * exp(-alpha1 * alpha2 / alpha_sum * (v1 - v2).squaredNorm());
}

int main() {
    ScalarType alpha1 = ScalarType(1.25);
    Vector<ScalarType, 3> v1{2, 5, -1};
    GaussBase<ScalarType> base1 = GaussBase<ScalarType>(v1, alpha1, 0, 0, 0);
    ScalarType alpha2 = ScalarType(0.76);
    Vector<ScalarType, 3> v2{-3, 6, 1};
    GaussBase<ScalarType> base2 = GaussBase<ScalarType>(v2, alpha2, 0, 0, 0);
    if (!scalarNear(base1.overlap(base2), overlap_1s_1s(alpha1, v1, alpha2, v2), 1E-14))
        return 1;
    return 0;
}
