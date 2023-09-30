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
#include <fstream>
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using namespace Physica::Utils;

const static char* data1 = "Structure\n"
                           "1.0\n"
                           "	   22.84502     0.00000     0.00000\n"
                           "	   -5.16077     4.54432     0.00000\n"
                           "	    0.00000     0.00000    12.00000\n"
                           " O H\n"
                           " 1 2\n"
                           "D\n"
                           "	    0.83817     0.05987     0.02619\n"
                           "	    0.97162     0.88296     0.02801\n"
                           "	    0.73317     0.52985     0.03865\n";

const static char* data2 = "Structure\n"
                           "1\n"
                           "-2.734285261           0 2.734285261\n"
                           "0            2.734285261 2.734285261\n"
                           "-2.734285261 2.734285261           0\n"
                           "Si\n"
                           "2\n"
                           "Direct\n"
                           "0.00 0.00 0.00\n"
                           "0.25 0.25 0.25";

Poscar<Scalar<Double>> readTest1() {
    auto tmp = TempFile("/tmp/tmpXXXXXX");
    std::ofstream os(tmp.getName());
    os << data1;
    os.close();

    Poscar<Scalar<Double>> poscar{};
    std::ifstream is(tmp.getName());
    is >> poscar;
    is.close();

    const auto& elements = poscar.getElementTypes();
    if (!(elements[0] == 8 && elements[1] == 1))
        exit(EXIT_FAILURE);
    const auto& numOfEachType = poscar.getNumOfEachType();
    if (!(numOfEachType[0] == 1 && numOfEachType[1] == 2))
        exit(EXIT_FAILURE);
    return poscar;
}

void readTest2() {
    auto tmp = TempFile("/tmp/tmpXXXXXX");
    std::ofstream os(tmp.getName());
    os << data2;
    os.close();

    Poscar<Scalar<Float>> poscar{};
    std::ifstream is(tmp.getName());
    is >> poscar;
    is.close();
}

int main() {
    {
        using PoscarType = Poscar<Scalar<Double>>;
        PoscarType poscar = readTest1();

        typename PoscarType::LatticeMatrix mat = poscar.getLattice();
        poscar.standrizeLattice();
        if (!matrixNear(mat, poscar.getLattice(), 1E-15))
            return 1;
    }
    {
        readTest2();
    }
    {
        using ScalarType = Scalar<Float>;
        using CrystalCellType = CrystalCell<ScalarType>;
        typename CrystalCellType::LatticeMatrix lattice{
            -4.6635062604325164,   -0.2499522611778955,    0.0000000000000000,
            -2.1629745970109657,   -4.1943944839773311,    0.0000000000000000,
            -0.2750800827878018,   -0.4169789280520980,   18.0000000000000000
        };
        typename CrystalCellType::PositionMatrix pos {
            0.4553508091084409,  0.3980437584135783,  0.1240303800896787,
            0.4937103263031835,  0.4030549988960055,  0.9488679230950712,
            0.5596918259357793,  0.8517822319914985,  0.1226285591691945,
            0.3686253245184842,  0.6403194088783717,  0.0163388989450929,
            0.5980496296529945,  0.8452761470000689,  0.9474667297556534,
            0.1076134753395919,  0.7454143691003363,  0.1235631952109221,
            0.6847635746074375,  0.9781822048654215,  0.0551553884196571,
            0.9457728898525665,  0.8447981726837335,  0.9479330783198919,
            0.6786065180112657,  0.9738018598142685,  0.1107906720462844,
            0.3747794285285615,  0.6414996922536187,  0.9607035572233092,
            0.7021261874138659,  0.9803177844871507,  0.9573706719037773,
            0.3512600170478342,  0.6392493670479714,  0.1141244566914832
        };
        const CrystalCellType unit({std::move(lattice), std::move(pos), CrystalCellType::Type::Direct}, {1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8});
        Poscar<ScalarType> poscar(unit);
        const auto& numOfEachType = poscar.getNumOfEachType();
        if (numOfEachType[0] != 8 || numOfEachType[1] != 4)
            return 1;
    }
    return 0;
}