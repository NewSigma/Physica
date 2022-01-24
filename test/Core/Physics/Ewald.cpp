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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/Ewald.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/ReciprocalCell.h"
#include "Physica/Core/Physics/PhyConst.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;

void VASPTest() {
    const double lengthInBohr = PhyConst<AU>::angstormToBohr(3);
    CrystalCell cell({lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {0.5, 0.5, 0.5}, {14});
    const auto repCell = cell.reciprocal();
    const auto ewald = Ewald<ScalarType>::energyIonIon(cell, repCell, {4});
    if (!scalarNear(ewald, ScalarType(PhyConst<AU>::eVToHartree(-108.95061336198556)), 1E-7))
        exit(EXIT_FAILURE);
}

/**
 * Reference:
 * [1] pyewald(github.com/lukeolson/pyewald)
 */
void madelungTest() {
    {
        const double lengthInBohr = PhyConst<AU>::angstormToBohr(5.6903014761756712);
        CrystalCell NaCl({lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {
                            0.0, 0.0, 0.0,
                            0.0, 0.5, 0.5,
                            0.5, 0.0, 0.5,
                            0.5, 0.5, 0.0,
                            0.5, 0.5, 0.5,
                            0.5, 0.0, 0.0,
                            0.0, 0.5, 0.0,
                            0.0, 0.0, 0.5
                        }, {1, 1, 1, 1, 1, 1, 1, 1});
        const auto ewald = Ewald<ScalarType>::energyIonIon(NaCl, NaCl.reciprocal(), {1, 1, 1, 1, -1, -1, -1, -1});
        const auto madelung = -(ewald / 4) * (lengthInBohr / 2); //We have 4x unit cell so energy is divided by 4
        if (!scalarNear(madelung, ScalarType(1.7475645946331822), 1E-7))
            exit(EXIT_FAILURE);
    }
    {
        const double lengthInBohr = 1;
        CrystalCell CsCl({lengthInBohr, 0, 0, 0, lengthInBohr, 0, 0, 0, lengthInBohr}, {
                            0.0, 0.0, 0.0,
                            0.5, 0.5, 0.5,
                        }, {1, 1});
        const auto ewald = Ewald<ScalarType>::energyIonIon(CsCl, CsCl.reciprocal(), {1, -1});
        const auto madelung = -ewald * (lengthInBohr * 0.5 * std::sqrt(3.0));
        if (!scalarNear(madelung, ScalarType(1.76267477307099), 1E-7))
            exit(EXIT_FAILURE);
    }
    {
        const double lengthInBohr = 0.5;
        CrystalCell ZnS({0, lengthInBohr, lengthInBohr, lengthInBohr, 0, lengthInBohr, lengthInBohr, lengthInBohr, 0}, {
                            0.0, 0.0, 0.0,
                            0.25, 0.25, 0.25,
                        }, {1, 1});
        const auto ewald = Ewald<ScalarType>::energyIonIon(ZnS, ZnS.reciprocal(), {1, -1});
        const auto madelung = -ewald * (lengthInBohr * 0.5 * std::sqrt(3.0));
        if (!scalarNear(madelung, ScalarType(1.63805505338879), 1E-8))
            exit(EXIT_FAILURE);
    }
}

int main() {
    VASPTest();
    madelungTest();
    return 0;
}
