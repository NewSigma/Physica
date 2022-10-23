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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/PWBaseWave.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;

int main() {
    const ScalarType cutOffE = ScalarType(0.5);
    PWBaseWave<ScalarType> wave(cutOffE, {1, 2, 3, 3, 2, 1, -1, 2, 6});
    ssize_t dimX = wave.getDimX();
    ssize_t dimY = wave.getDimY();
    ssize_t dimZ = wave.getDimZ();

    bool xWithinRange = wave.getKinetic(dimX, 0, 0) <= cutOffE;
    if (!xWithinRange)
        return 1;

    bool yWithinRange = wave.getKinetic(0, dimY, 0) <= cutOffE;
    if (!yWithinRange)
        return 1;

    bool zWithinRange = wave.getKinetic(0, 0, dimZ) <= cutOffE;
    if (!zWithinRange)
        return 1;

    bool xNotEnough = wave.getKinetic(dimX + 1, 0, 0) <= cutOffE;
    if (xNotEnough)
        return 1;

    bool yNotEnough = wave.getKinetic(0, dimY + 1, 0) <= cutOffE;
    if (yNotEnough)
        return 1;

    bool zNotEnough = wave.getKinetic(0, 0, dimZ + 1) <= cutOffE;
    if (zNotEnough)
        return 1;
    return 0;
}
