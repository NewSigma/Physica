/*
 * Copyright 2019 WeiBo He.
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
#include <iostream>
#include <Physica/Physica.h>
#include <Physica/Core/MultiPrecition/Scalar.h>

using namespace Physica::Core;

void printElements(const MultiScalar& n) {
    int size = n.getSize();
    for(int i = 0; i < size; ++i)
        std::cout << n[i] << ' ';
}

int main(int argc, char** argv) {
    initPhysica();
    /* test Pi */ {
        auto Pi = MathConst::getInstance().getPI();
        if(!(Pi.getSize() == 5
             && Pi[0] == 11424456171093639400UL
             && Pi[1] == 11820040416388919749UL
             && Pi[2] == 1376283091369227076
             && Pi[3] == 2611923443488327891
             && Pi[4] == 3)) {
            std::cout << "Current Pi: ";
            printElements(MathConst::getInstance().getPI());
            std::cout << '\n';
            return 1;
        }
    }
    /* test E */ {
        auto E = MathConst::getInstance().getE();
        if(!(E.getSize() == 4
             && E[0] == 7126689189968796226
             && E[1] == 13794904443024896967UL
             && E[2] == 13249961062380153450UL
             && E[3] == 2)) {
            std::cout << "Current E: ";
            printElements(MathConst::getInstance().getE());
            std::cout << '\n';
            return 1;
        }
    }
    deInitPhysica();
    return 0;
}