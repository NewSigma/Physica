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
#ifndef _Physica_C_ClimbMountainAlgorithm_H
#define _Physica_C_ClimbMountainAlgorithm_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    class HillClimbingAlgorithm {
    public:
        HillClimbingAlgorithm(MultiScalar* func(MultiScalar*), MultiScalar* x_initial, MultiScalar* stepSize);
        ~HillClimbingAlgorithm();
        void getExtremal();
        MultiScalar* getMinStep();
        void setMinStep(MultiScalar* minStep);
    private:
        MultiScalar* (*func)(MultiScalar*);
        MultiScalar* x_initial;
        MultiScalar* stepSize;
        MultiScalar* minStep;
    };
}

#endif