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
#include "Physica/Core/MultiPrecision/DiffTraceGuard.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"

using namespace Physica::Core;

using PlainScalar = Scalar<Double>;
using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse>;

int main() {
    SilveraGoldman<ScalarType> sg(1.0);
    {
        auto guard = DiffTraceGuard<PlainScalar>::make_guard();
        ScalarType r = 2.0;
        const ScalarType r2 = square(r);
        sg.pot_functor(0, 0, r, r2).reverse();
        const ScalarType f = -r.getTangent();
        const ScalarType f1 = sg.force_functor(0, 0, r, r2);
        if (!scalarNear(f, f1, 1E-15))
            return 1;
    }
    {
        auto guard = DiffTraceGuard<PlainScalar>::make_guard();
        ScalarType r = 2.0;
        const ScalarType r2 = square(r);
        sg.force_functor(0, 0, r, r2).reverse();
        const ScalarType fc = -r.getTangent();
        const ScalarType fc1 = sg.forceConst_functor(r, r2);
        std::cout << fc << ' ' << fc1 << std::endl;
        if (!scalarNear(fc, fc1, 1E-15))
            return 1;
    }
    return 0;
}
