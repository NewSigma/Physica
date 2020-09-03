/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_DIFFERENTIAL_H
#define PHYSICA_DIFFERENTIAL_H

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunction.h"

namespace Physica::Core {
    class Differential {
        TreeFunction<> func;
        MultiScalar at;
        MultiScalar stepSize;
    public:
        //Reference: Numerical Recipes in C++
        enum DifferentialMethod {
            DoublePoint,
            Forward,
            Backward
        };
        Differential(TreeFunction<> func, MultiScalar at
                , MultiScalar stepSize = MultiScalar(BasicConst::getInstance().getStepSize()));
        [[nodiscard]] MultiScalar solve(DifferentialMethod method);
    };
}

#endif
