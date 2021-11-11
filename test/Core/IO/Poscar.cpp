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
#include <unistd.h>
#include <fstream>
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Utils/TestHelper.h"

using namespace Physica::Core;
using namespace Physica::Utils;

const static char* data1 = "Structure	-97.8256	0.0003\n"
                           "1.0\n"
                           "	   22.84502     0.00000     0.00000\n"
                           "	   -5.16077     4.54432     0.00000\n"
                           "	    0.00000     0.00000    12.00000\n"
                           " 1 2\n"
                           "D\n"
                           "	    0.83817     0.05987     0.02619\n"
                           "	    0.97162     0.88296     0.02801\n"
                           "	    0.73317     0.52985     0.03865\n";

Poscar readTest() {
    const char* tempFile = tmpnam(nullptr);
    std::ofstream os(tempFile);
    os << data1;
    os.close();

    Poscar poscar{};
    std::ifstream is(tempFile);
    is >> poscar;
    is.close();

    unlink(tempFile);

    const auto& numOfEachType = poscar.getNumOfEachType();
    if (!(numOfEachType[0] == 1 && numOfEachType[1] == 2))
        exit(EXIT_FAILURE);
    return poscar;
}

int main() {
    Poscar poscar = readTest();

    typename Poscar::LatticeMatrix mat = poscar.getLattice();
    poscar.standrizeLattice();
    if (!matrixNear(mat, poscar.getLattice(), 1E-15))
        return 1;
    return 0;
}