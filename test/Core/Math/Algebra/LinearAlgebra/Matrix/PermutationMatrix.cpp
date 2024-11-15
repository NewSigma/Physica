/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermutationMatrix.h"
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;
using ScalarType = float64;

int main() {
    const PermutationMatrix<ScalarType> perm({0, 3, 2, 1});
    const VectorND<ScalarType> v{0, 1, 2, 3};
    const VectorND<ScalarType> perm_v = perm * v;
    {
        const VectorND<ScalarType> answer{0, 3, 2, 1};
        if (!vectorNear(perm_v, answer, 1E-15))
            return 1;
    }
    const auto inv = perm.inverse();
    const VectorND<ScalarType> v1 = inv * perm_v;
    if (!vectorNear(v, v1, 1E-15))
        return 1;
    return 0;
}
