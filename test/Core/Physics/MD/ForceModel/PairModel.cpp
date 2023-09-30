/*
 * Copyright 2023 WeiBo He.
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
#include "Physica/Core/MultiPrecision/Differentiable.h"
#include "Physica/Core/Physics/MD/ForceModel/LJModel.h"

using namespace Physica::Core;

using ScalarType = Differentiable<Scalar<Double>, DiffMode::Reverse>;
using PosScalarType = ScalarType;

int main() {
    LJModel<ScalarType, PosScalarType> lj(1.0, 1.0);
    ScalarType r = 1.0;
    lj.pot_functor(0, 0, r, square(r)).reverse();
    const ScalarType f = r.getTangent();
    const ScalarType f1 = lj.force_functor(0, 0, r, square(r));
    if (!scalarNear(f, -f1, 1E-15))
        return 1;
    return 0;
}
