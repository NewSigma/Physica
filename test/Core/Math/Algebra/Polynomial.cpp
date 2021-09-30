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
#include "Physica/Core/Math/Algebra/BasicAlgebra/Polynomial.h"

using namespace Physica::Core;

template<class VectorType>
bool testPolyRoot(const RValueVector<VectorType>& coeffs, double precision) {
    using ScalarType = typename VectorType::ScalarType;
    const Polynomial<ScalarType, VectorType::SizeAtCompile> poly(coeffs);
    auto roots = polyRoot(poly);
    for (const auto& root : roots) {
        auto result = poly(root);
        if (!scalarNear(result, decltype(result)::Zero(), precision))
            return false;
    }
    return true;
}

int main() {
    using ScalarType = Scalar<Double, false>;
    Vector<ScalarType, 6> coeffs{1, 2, 3, 4, 5, 6};
    if (!testPolyRoot(coeffs, 1E-11))
        return 1;
    return 0;
}
